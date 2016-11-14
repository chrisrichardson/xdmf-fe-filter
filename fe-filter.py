## execfile("/home/chris/code/FEniCS/src/kitware/fe-filter.py")

def make_outcell_from_tetrahedron(in_vals, element_type,
                                  value_indices, in_cell,
                                  ds_out, out_vals):

    cells_points = in_cell.GetPoints()

    index_nc = len(value_indices)
    # Work out number of components (scalar/vector/tensor)
    ncomp = 1
    if element_type[:2] == "CG" or element_type[:2] == "DG":
        d = int(element_type[2:])
        d = (d + 1)*(d + 2)*(d + 3)/6
        assert index_nc%d == 0
        ncomp = index_nc/d

    out_vals.SetNumberOfComponents(ncomp)

    if element_type == "DG1" or element_type == "CG1":
        # CG1 should be simplified (no need to rebuild mesh)
        ptlist = [None, None, None, None]
        scnt = ds_out.GetPoints().GetNumberOfPoints()
        for ix in range(0, 4):
            coord = cells_points.GetPoint(ix)
            ds_out.GetPoints().InsertNextPoint(coord)

            # Scalar (one value), Vector (two or three values)
            for i in range(ncomp):
                index = value_indices[ix + i*d]
                v = in_vals.GetTuple(index)[0]
                out_vals.InsertNextValue(v)

            ptlist[ix] = scnt
            scnt = scnt+1
        ds_out.InsertNextCell(vtk.VTK_TETRA, 4, ptlist)

    elif element_type == "CG2" or element_type == "DG2":
        # Represent CG2/DG2 on a quadratic mesh
        ptlist = [None, None, None, None, None, None, None, None, None, None]
        scnt = ds_out.GetPoints().GetNumberOfPoints()
        for ix in range(0, 10):
            if ix < 4:
                c = cells_points.GetPoint(ix)
            else:
                c0 = cells_points.GetPoint((ix - 1)%3)
                iy = ix%3
                if ix > 6:
                    iy = 3
                c1 = cells_points.GetPoint(iy)
                c = [0.5*(c0[0] + c1[0]),
                     0.5*(c0[1] + c1[1]),
                     0.5*(c0[2] + c1[2])]
            ds_out.GetPoints().InsertNextPoint(c)

            qmap = [0, 1, 2, 3, 9, 6, 8, 7, 5, 4]
            for i in range(ncomp):
                index = value_indices[qmap[ix] + i*d]
                v = in_vals.GetTuple(index)[0]
                out_vals.InsertNextValue(v)

            ptlist[ix] = scnt
            scnt = scnt+1
        ds_out.InsertNextCell(vtk.VTK_QUADRATIC_TETRA, 10, ptlist)
    else:
        print "Cannot yet represent element type: ", element_type
    return

def make_outcell_from_triangle(in_vals, element_type,
                               value_indices, in_cell,
                               ds_out, out_vals):

    index_nc = len(value_indices)
    # Work out number of components (scalar/vector/tensor)
    ncomp = 1
    if element_type[:2] == "CG" or element_type[:2] == "DG":
        d = int(element_type[2:])
        d = (d + 1)*(d + 2)/2
        assert index_nc%d == 0
        ncomp = index_nc/d

    if (element_type == "RT1"):
        ncomp = 3

    out_vals.SetNumberOfComponents(ncomp)

    cells_points = in_cell.GetPoints()

    if element_type == "DG1" or element_type == "CG1":
        # CG1 should be simplified (no need to rebuild mesh)
        ptlist = [None, None, None]
        scnt = ds_out.GetPoints().GetNumberOfPoints()
        for ix in range(0, 3):
            coord = cells_points.GetPoint(ix)
            ds_out.GetPoints().InsertNextPoint(coord)

            # Scalar (one value), Vector (two or three values)
            for i in range(ncomp):
                index = value_indices[ix + i*3]
                v = in_vals.GetTuple(index)[0]
                out_vals.InsertNextValue(v)

            ptlist[ix] = scnt
            scnt = scnt+1
        ds_out.InsertNextCell(vtk.VTK_TRIANGLE, 3, ptlist)

    elif element_type == "CG2" or element_type == "DG2":
        # Represent CG2/DG2 on a quadratic mesh
        ptlist = [None, None, None, None, None, None]
        scnt = ds_out.GetPoints().GetNumberOfPoints()
        for ix in range(0, 6):
            if ix < 3:
                c = cells_points.GetPoint(ix)
            else:
                c0 = cells_points.GetPoint(ix - 3)
                c1 = cells_points.GetPoint((ix - 2)%3)
                c = [0.5*(c0[0] + c1[0]),
                     0.5*(c0[1] + c1[1]),
                     0.5*(c0[2] + c1[2])]
            ds_out.GetPoints().InsertNextPoint(c)

            qmap = [0, 1, 2, 5, 3, 4]
            for i in range(ncomp):
                index = value_indices[qmap[ix] + i*6]
                v = in_vals.GetTuple(index)[0]
                out_vals.InsertNextValue(v)

            ptlist[ix] = scnt
            scnt = scnt+1
        ds_out.InsertNextCell(vtk.VTK_QUADRATIC_TRIANGLE, 6, ptlist)

    elif element_type[:2] == "CG":
        # Downsample CG* to CG1 by taking the first three values...
        # Could be much improved, e.g. by creating more cells...
        ptlist = [None, None, None]
        scnt = ds_out.GetPoints().GetNumberOfPoints()
        for ix in range(0, 3):
            coord = cells_points.GetPoint(ix)
            ds_out.GetPoints().InsertNextPoint(coord)

            for i in range(nv):
                index = value_indices[ix + i*d]
                v = in_vals.GetTuple(index)[0]
                out_vals.InsertNextValue(v)

            ptlist[ix] = scnt
            scnt = scnt+1
        ds_out.InsertNextCell(vtk.VTK_TRIANGLE, 3, ptlist)

    elif element_type == "RT1":

        assert index_nc == 3
        ptlist = [None, None, None]
        scnt = ds_out.GetPoints().GetNumberOfPoints()

        point_ids = [in_cell.GetPointIds().GetId(i) for i in range(3)]
        coord = [cells_points.GetPoint(i) for i in range(3)]
        print point_ids

        # Calculate normals
        n = [None, None, None]
        for ix in range(3):
            ix1 = (ix+1)%3
            n[ix] = [coord[ix1][1] - coord[ix][1],
                     coord[ix][0] - coord[ix1][0]]
            n[ix] /= sqrt(n[ix][0]**2 + n[ix][1]**2)
            if (point_ids[ix] > point_ids[ix1]):
                n[ix] *= -1

        for ix in range(0, 3):
            pt = [0.5*(coord[ix][0] + coord[(ix+1)%3][0]),
                  0.5*(coord[ix][1] + coord[(ix+1)%3][1]),
                  0.0]
            ds_out.GetPoints().InsertNextPoint(pt)

            index = value_indices[(ix+2)%3]
            v = in_vals.GetTuple(index)[0]
            out_vals.InsertNextValue(v*n[ix][0])
            out_vals.InsertNextValue(v*n[ix][1])
            out_vals.InsertNextValue(0.0)

            ptlist[ix] = scnt
            scnt = scnt+1
        ds_out.InsertNextCell(vtk.VTK_POLY_VERTEX, 3, ptlist)
    else:
        print "Cannot yet represent element type: ", element_type
    return

def make_outcell_from_incell(in_vals, in_cell, element_type,
                             value_indices, ds_out, out_vals):

    global make_outcell_from_triangle
    global make_outcell_from_tetrahedron

    in_cell_type = in_cell.GetCellType()

    if in_cell_type == vtk.VTK_TRIANGLE:
        make_outcell_from_triangle(in_vals, element_type, value_indices, in_cell, ds_out, out_vals)
    elif in_cell_type == vtk.VTK_TETRA:
        make_outcell_from_tetrahedron(in_vals, element_type, value_indices, in_cell, ds_out, out_vals)
    else:
        print "Cannot yet represent cell type: ", in_cell_type

def traverse(ds_in, ds_out):
    global make_outcell_from_incell

    # Get a hold of arrays I am going to transcribe
    # Should be one set of CellData (indices)
    # and one set of FieldData (values)
    assert ds_in.GetCellData().GetNumberOfArrays() == 1
    in_idx_name = ds_in.GetCellData().GetArrayName(0)
    index_array = ds_in.GetCellData().GetArray(in_idx_name)
    index_nc = index_array.GetNumberOfComponents()

    assert ds_in.GetFieldData().GetNumberOfArrays() == 1
    in_vals_name = ds_in.GetFieldData().GetArrayName(0)
    in_vals = ds_in.GetFieldData().GetArray(in_vals_name)

    # Get the element type and check it matches
    assert in_idx_name.split("_")[-1] == "idx"
    element_type = in_idx_name.split("_")[-2]
    assert in_vals_name.split("_")[-1] == "val"
    assert in_vals_name.split("_")[-2] == element_type

    print "Element Type = ", element_type

    out_vals = in_vals.NewInstance()
    out_vals.SetName(in_vals.GetName())

    # For DG0 elements, just add a Cell Array and return
    if element_type == "DG0":
        out_vals.SetNumberOfComponents(index_nc)
        out_vals.Allocate(ds_in.GetNumberOfCells())
        ds_out.GetCellData().AddArray(out_vals)
        for x in range(ds_in.GetNumberOfCells()):
            for index in index_array.GetTuple(x):
                v = in_vals.GetTuple(int(index))[0]
                out_vals.InsertNextValue(v)
        # Remove FieldData - no longer needed
        ds_out.GetFieldData().RemoveArray(in_vals_name)
        return
    else:
        ds_out.Allocate(ds_in.GetNumberOfCells())
        ds_out.SetPoints(vtk.vtkPoints())
        ds_out.GetPointData().AddArray(out_vals)

    for x in range(0, ds_in.GetNumberOfCells()):
        in_cell = ds_in.GetCell(x)

        # Indices to quantities for this cell
        value_indices = [int(i) for i in index_array.GetTuple(x)]

        # Now make up new cells for that input cell
        make_outcell_from_incell(in_vals, in_cell, element_type,
                                 value_indices, ds_out, out_vals)

    # Remove FieldData - no longer needed
    ds_out.GetFieldData().RemoveArray(in_vals_name)

ds_in = self.GetInput()
ds_out = self.GetOutput()
traverse(ds_in, ds_out)
