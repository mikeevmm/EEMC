#!/usr/bin/env python3

import matplotlib.pyplot as plt
import os 

for file in ('/ex2-heh', 'ex2-h'):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    infile = dir_path + file + '.out

    pts = []
    min_val = None

    with open(infile) as infile:
        for line in infile.read().splitlines():
            if line.startswith('#'):
                continue
            x,y = tuple(map(lambda x: float(x), filter(lambda x: x, line.split())))
            pts.append((x,y))

            if min_val is None or y < min_val[1]:
                min_val = (x,y)

    x_pts = tuple(map(lambda x: x[0], pts))
    y_pts = tuple(map(lambda x: x[1], pts))
    plt.plot(x_pts, y_pts, '+', label='Calculated points')

    plt.plot(min_val[0], min_val[1], 'o', color='r', label='Min@({:.4f},{:.4f})'.format(min_val[0], min_val[1]))

    plt.xlabel('Distance (Angstrom)')
    plt.ylabel('Energy (Hartree)')
    plt.title('HeH+ E(d)')
    plt.legend(loc='best')

    #plt.savefig(file + '.png')

    plt.show()
