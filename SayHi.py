#!/usr/bin/env python3

import socket
import threading
import time
import random

class Channel(object):
	WHOAMI = 0
	PROBE = 1

	def __init__(self, in_port):
		self.channel = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
		self.channel.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
		self.channel.bind(("", in_port))
		self.channel.settimeout(0)
		self.channel_lock = threading.Lock()
		self.port = in_port

		self.addr_secret = random.getrandbits(32).to_bytes(4, byteorder='big')
		self.addr = None

	def __enter__(self):
		self.listen_thread = threading.Thread(target=self._thread_listen_)
		self.listen_alive = True
		self.listen_thread.start()
		self._find_addr_()
		return self

	def __exit__(self, exc_type, exc_value, traceback):
		self.channel.close()
		self.listen_alive = False
		self.listen_thread.join()

	def _thread_listen_(self):
		while (self.listen_alive):
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
				self._on_msg_(msg, addr)
			time.sleep(0.1)
	
	def _on_msg_(self, msg, addr):
		msg = bytearray(msg)
		id = int.from_bytes(msg[:2], byteorder='big')
		payload = msg[2:]

		if id == self.WHOAMI:
			if payload == self.addr_secret:
				print('Found self addr:', addr)
				self.addr = addr		

	def _find_addr_(self):
		id = self.identifier(self.WHOAMI) 	# 2 bytes
		random_key = self.addr_secret		# 4 bytes
		byte_digest = id + random_key
		self.broadcast(byte_digest)
	
	def identifier(self, id_int):
		return id_int.to_bytes(2, byteorder='big')

	def broadcast(self, msg):
		self.send(msg, '<broadcast>', self.port)

	def send(self, msg, addr, addr_port):
		self.channel_lock.acquire()
		self.channel.sendto(msg, (addr, addr_port))
		self.channel_lock.release()


with Channel(6969) as channel:
	while True:
		pass