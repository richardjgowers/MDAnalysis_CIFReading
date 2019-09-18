"""Example of converting a cif file to PDB using MDAnalysis

"""
import MDAnalysis as mda
import mdaCIF

u = mda.Universe('IRMOF-1.cif')

u.atoms.write('out.pdb')

