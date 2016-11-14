
# xdmf-fe-filter

ParaView filter to read Finite Element data in XDMF format

## Format
At present, the XDMF contains an Unstructured Grid, i.e. mesh topology and geometry, and two *Attributes*, which
must be named as *xxx_ELEMENT_idx* and *xxx_ELEMENT_val*, where *xxx* is user defined, and *ELEMENT* is one of:

 * CG1, CG2
 * DG0, DG1, DG2
 * RT1

The *xxx_ELEMENT_val* Attribute is a 1D float array, and the *xxx_ELEMENT_idx* is a 2D UInt array, with index values
for the nodal values in each cell. Hence the first dimension must be the same as the number of cells.
In order to read correctly, use the Xdmf3 reader in ParaView.

A script `generate.py` is provided, which uses experimental XDMF output in FEniCS to generate the XDMF files.

## Testing in ParaView

 1. Open the XDMF file with the Xdmf3 reader
 2. Add a filter "ProgrammableFilter" and in the box paste: `execfile("fe-filter.py")` - you may need to add the full path.
 3. Select the output field of the filter to display and process as usual

## Future Work

Add support for more elements (higher order, BDM, etc) , and more cell types (e.g. 1D, quadratic meshes)
