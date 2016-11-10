#!/usr/bin/env python
#
# Code for calculating angles against BioPython PDB module.
# Also does some formatting.
#
# When called as a script dump out the angles for a specified PDB
# file and chain.

from Bio.PDB.Vector import calc_dihedral, calc_angle
from math import pi

def deg(theta): return theta*180.0/pi

def v_(pp, i, kind):
    """shorthand to extract vector"""
    return pp[i][kind].get_vector()

def phi(pp, i):
    """phi angle for peptide index i"""
    cp=v_(pp, i-1, 'C')
    n=v_(pp, i, 'N')
    ca=v_(pp, i, 'CA')
    c=v_(pp, i, 'C')
    return deg(calc_dihedral(cp, n, ca, c))

def psi(pp, i):
    """psi angle for peptide index i"""
    n=v_(pp, i, 'N')
    ca=v_(pp, i, 'CA')
    c=v_(pp, i, 'C')
    nn=v_(pp, i+1, 'N')
    return deg(calc_dihedral(n, ca, c, nn))

def omega(pp, i):
    """omega angle for peptide index i"""
    ca=v_(pp, i, 'CA')
    c=v_(pp, i, 'C')
    nn=v_(pp, i+1, 'N')
    can=v_(pp, i+1, 'CA')
    return deg(calc_dihedral(ca, c, nn, can))

def ca_angle(pp, i):
    """alpha carbon bond angle for peptide index i"""
    n=v_(pp, i, 'N')
    ca=v_(pp, i, 'CA')
    c=v_(pp, i, 'C')
    return deg(calc_angle(n,ca,c))

def c_angle(pp, i):
    """carbon bond angle for peptide index i"""
    ca=v_(pp, i, 'CA')
    c=v_(pp, i, 'C')
    nn=v_(pp, i+1, 'N')
    return deg(calc_angle(ca,c,nn))

def n_angle(pp, i):
    """nitrogen bond angle for peptide index i"""
    cp=v_(pp, i-1, 'C')
    n=v_(pp, i, 'N')
    ca=v_(pp, i, 'CA')
    return deg(calc_angle(cp, n, ca))

fns=(phi, psi, omega, n_angle, ca_angle, c_angle )
def get_angle_table(pp, i, lng=9):
    """length lng table of angles starting at peptide i"""
    tab=[]
    for i in range(idx, lng):
        tab.append(fn(pp,i) for fn in fns)
    return tuple(tab)

def p(n): print ("%.1f" % n).rjust(6),
def q(n): print " " + "-"*(n-1),

def print_peptide(pp, i):
    lng=len(pp)
    # Phi
    if i>0:
        p(phi(pp, i))
    else:
        q(6)
    # Psi and Omega
    if i<lng-1:
        p(psi(pp, i))
        p(omega(pp, i))
    else:
        q(6); q(6)
    # N
    if i>0:
        p(n_angle(pp, i))
    else:
        q(6)
    # CA
    p(ca_angle(pp, i))
    # C
    if i<lng-1:
        p(c_angle(pp, i))
    else:
        q(6)
    print ""

header="#  idx  resid   name   phi    psi  omega      n      c     ca"
   
def print_pp(pp, offset=0, lng=None):
    """print a polypeptide pp, perhaps with offset offset and length lng"""
    if not lng:
        lng=len(pp)    
    #print header
    for i in range(offset, offset+lng):
        # offset, sequence ID and name
        print ("%d" % (int(i)-int(offset))).rjust(6),
        # print ("%d" % pp[i].get_id()[1]).rjust(6),
        # print (pp[i].get_resname()).rjust(6),
        print_peptide(pp, i)

if __name__ == "__main__":
    from sys import argv
    from pdb_extract import get_pps
    
    pps=get_pps(argv[1], argv[2])
    for pp in pps:
        print header
        for i in range(0, len(pp)):
            # offset, sequence ID and name
            print ("%d" % i).rjust(6),
            print ("%d" % pp[i].get_id()[1]).rjust(6),
            print (pp[i].get_resname()).rjust(6),            
            print_peptide(pp, i)


