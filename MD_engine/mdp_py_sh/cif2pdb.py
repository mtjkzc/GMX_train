### CIF to PDB converter ###
# Author: Matej Kožić | mkozic@chem.pmf.hr
# End of comment

filename = "cif_filename"
savename = "output_pdb_filename"

import Bio.PDB as bp

parser = bp.MMCIFParser()
cif_structure = parser.get_structure(" ", filename)

io = bp.PDBIO()
io.set_structure(cif_structure)
io.save(savename)

# END
