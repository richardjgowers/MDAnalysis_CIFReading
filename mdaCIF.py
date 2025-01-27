import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.base import SingleFrameReaderBase
from MDAnalysis.topology.base import TopologyReaderBase
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    Atomnames,
    Charges,
    Masses,
    Resids,
    Resnums,
    Segids,
)

from openbabel import pybel
import openbabel
fillUC = openbabel.OBOp.FindType('fillUC')

def unpack_symmetry(mol):
    mol.unitcell.FillUnitCell(mol.OBMol)


class CIFReader(SingleFrameReaderBase):
    """Read coordinate information from a CIF file

    Reads
     - coordinates
     - unitcell information
    """
    format = 'CIF'

    def _read_first_frame(self):
        pbmol = next(pybel.readfile('cif', self.filename))
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

        ts.dimensions = A, B, C, alpha, beta, gamma

        return ts


class CIFParser(TopologyReaderBase):
    """Parses a CIF file and creates a MDAnalysis topology

    Reads
     - name
     - mass
     - charge
    """
    format = 'CIF'

    def parse(self, **kwargs):
        pbmol = next(pybel.readfile('cif', self.filename))
        unpack_symmetry(pbmol)
        
        natoms = len(pbmol.atoms)
        
        names = np.zeros(natoms, dtype=object)
        charges = np.zeros(natoms, dtype=np.float32)
        masses = np.zeros(natoms, dtype=np.float64)

        for i, atom in enumerate(pbmol.atoms):
            names[i] = atom.type
            masses[i] = atom.exactmass
            charges[i] = atom.partialcharge

        attrs = [
            Atomnames(names),
            Charges(charges),
            Masses(masses),
            Resids(np.array([1])),
            Resnums(np.array([1])),
            Segids(np.array(['SYSTEM'], dtype=object)),
        ]

        return Topology(natoms, 1, 1, attrs=attrs)
