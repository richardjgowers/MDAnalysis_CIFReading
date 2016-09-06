import MDAnalysis as mda
from MDAnalysis.coordinates.base import SingleFrameReader
from MDAnalysis.topology.base import TopologyReader
from MDAnalysis.core.AtomGroup import Atom

import pybel
import openbabel
fillUC = openbabel.OBOp.FindType('fillUC')

def unpack_symmetry(mol):
    mol.unitcell.FillUnitCell(mol.OBMol)


class CIFReader(SingleFrameReader):
    format = 'CIF'
    
    def _read_first_frame(self):
        pbmol = pybel.readfile('cif', self.filename).next()
        unpack_symmetry(pbmol)
        
        self.n_atoms = n_atoms = len(pbmol.atoms)
        self.ts = ts = self._Timestep(n_atoms, **self._ts_kwargs)
        for i, at in enumerate(pbmol.atoms):
            ts.positions[i] = at.coords
            
        A = pbmol.unitcell.GetA()
        B = pbmol.unitcell.GetB()
        C = pbmol.unitcell.GetC()
        alpha = pbmol.unitcell.GetAlpha()
        beta = pbmol.unitcell.GetBeta()
        gamma = pbmol.unitcell.GetGamma()
        
        ts._unitcell[:] = A, B, C, alpha, beta, gamma
            
        return ts


class CIFParser(TopologyReader):
    format = 'CIF'
    
    def parse(self):
        pbmol = pybel.readfile('cif', self.filename).next()
        unpack_symmetry(pbmol)
        
        natoms = len(pbmol.atoms)
        
        atoms = []
        for i, atom in enumerate(pbmol.atoms):
            name = elem = atom.type
            mass = atom.exactmass
            charge = atom.partialcharge
            atoms.append(Atom(i, name, elem,
                              'Res1', 1, 'SegA',
                              mass, charge, universe=self._u))
        return {'atoms': atoms}
