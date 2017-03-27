# Reading cif files with MDAnalysis

Reads cif files and unpacks all symmetries

## Requirements
 * [MDAnalysis](https://anaconda.org/MDAnalysis/mdanalysis)
 * [pybel](https://anaconda.org/patrickfuller/openbabel)

## Usage

Create a MDAnalysis Universe.  For how to use this, [see this tutorial](http://www.mdanalysis.org/MDAnalysisTutorial/)

```python
import MDAnalysis as mda
from mdaCIF import CIFReader, CIFParser

u = mda.Universe('yourcif.cif')
```

Convert cif into another format:

```python
import MDAnalysis as mda
from mdaCIF import CIFReader, CIFParser

u = mda.Universe('yourcif.cif')

# Can also write as xyz, gro etc
u.atoms.write('newfile.pdb')
```
