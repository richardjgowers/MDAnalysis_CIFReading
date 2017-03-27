"""Example of converting a cif file to PDB using MDAnalysis

"""
from mdaCIF import CIFReader, CIFParser
import MDAnalysis as mda

u = mda.Universe('IRMOF-1.cif')

u.atoms.write('out.pdb')

