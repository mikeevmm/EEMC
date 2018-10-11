import socket
import threading
import time
import random
import re

class Channel(object):
	WHOAMI = 0
	PROBEOUT = 1
	PROBEIN = 2
	TEXTMSG = 3

	"""
	Payload structures:
	WHOAMI:
		[2 bytes; payload ID] [4 bytes; random identifier] 
		Broadcast message to determine IP in subnet

	PROBEOUT:
		[2 bytes; payload ID] [1 byte; broadcast mode --- should include self in PROBEIN] [2 bytes (unsigned); port listening for PROBEIN]
		Request for the contact list of the peer
	
	PROBEIN:
		[2 bytes; payload ID] [1 bytes (unsigned); uint count of incoming peer address count ]
			[remainder; (ip address, port, human_name) in groups of 4+2+16 bytes]
		Response to PROBEOUT request from a queried peer
	
	TEXTMSG:
		[2 bytes; payload ID] [1 bytes (unsigned); private or public (broadcast) message] [remainder; text message]
	"""

	def __init__(self, in_port, human_name):
		self.channel = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
		self.channel.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
		self.channel.bind(("", in_port))
		self.channel.settimeout(0)
		self.channel_lock = threading.Lock()
		self.port = in_port
		self.human_name = human_name
		if len(human_name) > 8:
			print("Human name truncated to 8 characters")
			self.human_name = self.human_name[:8]

		self.addr_secret = random.getrandbits(32).to_bytes(4, byteorder='big')
		self.addr = None
		self.addr_human = None

		self.contact_list = set()
		self.human_to_contact = dict()

	def __enter__(self):
		self.alive = True
		
		self.listen_thread = threading.Thread(target=self._thread_listen_)
		self.listen_thread.start()

		while self.addr is None:
			id = self._get_identifier_bytes_(self.WHOAMI) 	# 2 bytes
			random_key = self.addr_secret		# 4 bytes
			byte_digest = id + random_key
			self._broadcast_(byte_digest)
			time.sleep(1)

		self.probe_thread = threading.Thread(target=self._thread_probe_)
		self.probe_thread.start()

		return self

	def __exit__(self, exc_type, exc_value, traceback):
		self.alive = False
		self.listen_thread.join()
		self.probe_thread.join()
		self.channel.close()

	def _thread_listen_(self):
		while (self.alive):
			self.channel_lock.acquire()
			msg, addr = None, None
			try:
				msg, addr = self.channel.recvfrom(1024)
				addr, port = addr
			except:
				pass # Assume EWOULDBLOCK
			finally:
				self.channel_lock.release()
			if None not in (msg, addr):
				self._on_msg_(msg, addr, port)
	
	def _thread_probe_(self):
		while(self.alive):
			if len(self.contact_list) == 0:
				print("No contacts, bootstrapping through <broadcast>")
				self._send_probeout_('', self.port, True)
			time.sleep(3)

	def _on_msg_(self, msg, addr, port):
		msg = bytearray(msg)
		id = int.from_bytes(msg[:2], byteorder='big', signed = False)
		payload = msg[2:]

		if id == self.WHOAMI: # Broadcast response to WHOAMI query
			if payload == self.addr_secret:
				print('Found self addr:', addr)
				self.addr_human = addr
				self.addr = bytes(map(int, addr.split('.')))

		elif id == self.PROBEOUT: # External request for contacts
			if addr == self.addr_human:
				return
			print("Recieved PROBEOUT")
			broadcast_mode = payload[0]
			if not broadcast_mode and len(self.contact_list) == 0:
				return
			return_list = self.contact_list
			if broadcast_mode:
				return_list.add((self.addr, self.port, self.human_name))
			if addr in return_list:
				return_list.remove(addr)
			
			requesting_port = int.from_bytes(payload[1:3], byteorder='big', signed=False)
			contacts_len_digest = int(len(return_list)).to_bytes(1, byteorder='big', signed=False)
			contacts_digest = bytearray(*((x[0] + int(x[1]).to_bytes(2, byteorder='big', signed=False) + bytes(self.human_name, 'utf-8')) for x in return_list))
			response_id = self._get_identifier_bytes_(self.PROBEIN)
			
			print("Got PROBEOUT request from {}, sending PROBEIN back".format(addr))
			full_response = response_id + contacts_len_digest + contacts_digest
			self._send_(full_response, addr, requesting_port) # Send PROBEIN back

		elif id == self.PROBEIN: # Response to contacts request
			print("Recieved PROBEIN")
			contact_count = payload[0]

			contacts = set()
			for i in range(1, len(payload[1:]), 4+2+8):
				contact_ip = bytes(payload[i:i+4])
				contact_port = int.from_bytes(payload[i+4:i+6], byteorder='big', signed=False)
				contact_name = bytes(payload[i+6:i+22]).decode('utf-8')
				contacts.add((contact_ip, contact_port, contact_name))
				
			new_contacts = contacts - self.contact_list
			self.contact_list = self.contact_list.union(new_contacts)

			print("Got PROBEIN response with {} new elements.".format(contact_count))
			for contact in new_contacts:
				contact_addr, contact_port, contact_name = contact
				# Save human reference
				print("Got new contact: (%s, %s, %s)" % (contact_addr, contact_port, contact_name))
				self.human_to_contact[contact_name] = (contact_addr, contact_port)
				# Send PROBEOUT to new contacts
				self._send_probeout_(contact_addr, contact_port)

		elif id == self.TEXTMSG:
			is_public_msg = payload[0]
			msg_content = payload[1:].decode('utf-8')
			print(msg_content)

	def _send_probeout_(self, to_addr, to_port, broadcast=False):
		probeout_id = self._get_identifier_bytes_(self.PROBEOUT)
		broadcast_mode = b'\x01' if broadcast else b'\x00'
		listening_on = self.port.to_bytes(2, byteorder='big', signed=False)
		msg_digest = probeout_id + broadcast_mode + listening_on
		if broadcast:
			self._broadcast_(msg_digest, to_port)
		else:
			self._send_(msg_digest, to_addr, to_port)
	
	def _get_identifier_bytes_(self, id_int):
		return id_int.to_bytes(2, byteorder='big', signed = False)

	def _broadcast_(self, msg, port = None):
		if port is None:
			port = self.port
		self._send_(msg, '<broadcast>', port)

	def _send_(self, msg, addr, addr_port):
		if type(addr) in (bytes, bytearray):
			addr = '.'.join(str(x) for x in addr)
		if type(addr_port) in (bytes, bytearray):
			addr_port = int.from_bytes(addr_port, byteorder='big', signed=False)
		self.channel_lock.acquire()
		self.channel.sendto(msg, (addr, addr_port))
		self.channel_lock.release()
	
	def broadcast_text(self, msg):
		if len(msg) == 0:
			return
		textmsg_id = self._get_identifier_bytes_(self.TEXTMSG)
		is_broadcast = b'\x01'
		text_bytes = str(msg).encode('utf-8')
		self._broadcast_(textmsg_id + is_broadcast + text_bytes)

	def send_text_to_human(self, msg, human):
		if len(msg) == 0:
			return
		if human not in self.human_to_contact:
			print("Unknown %s" % human)
			return
		addr, port = self.human_to_contact[human]
		textmsg_id = self._get_identifier_bytes_(self.TEXTMSG)
		is_broadcast = b'\x00'
		text_bytes = str(msg).encode('utf-8')
		self._send_(text_bytes, addr, port)

mention_reg = re.compile(r"[^\\]?@([\w_-]+)")

name = input("NAME>").strip().lower()
while mention_reg.match('@' + name) is None:
	print('Invalid username, use letters, numbers, - and _.')
	name = input("NAME>").strip().lower()

with Channel(1337, name) as channel:
	while True:
		msg = input(">").strip()
		sendto = tuple(map(lambda x: x.group(1), mention_reg.finditer(msg)))
		if len(sendto) == 0:
			channel.broadcast_text('%s (Public)> %s' % (name, msg))
		else:
			for dest in sendto:
				channel.send_text_to_human('%s (Private)> %s' % (name, msg), dest)