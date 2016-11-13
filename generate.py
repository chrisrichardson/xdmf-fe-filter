"""
Script to generate a range of XDMF files with Finite Element "Index + Data" format
"""

from dolfin import *

def output_scalar_xdmf(mesh, scalar_expression, file_suffix):
    Q = FunctionSpace(mesh, "DG", 0)
    F = Function(Q)
    F.interpolate(Expression(scalar_expression, degree=1))
    xdmf = XDMFFile("scalar_dg0_" + file_suffix + ".xdmf")
    xdmf.write_experimental(F, XDMFFile.Encoding_ASCII)

    Q = FunctionSpace(mesh, "DG", 1)
    F = Function(Q)
    F.interpolate(Expression(scalar_expression, degree=1))
    xdmf = XDMFFile("scalar_dg1_" + file_suffix + ".xdmf")
    xdmf.write_experimental(F, XDMFFile.Encoding_ASCII)

    Q = FunctionSpace(mesh, "DG", 2)
    F = Function(Q)
    F.interpolate(Expression(scalar_expression, degree=1))
    xdmf = XDMFFile("scalar_dg2_" + file_suffix + ".xdmf")
    xdmf.write_experimental(F, XDMFFile.Encoding_ASCII)

    Q = FunctionSpace(mesh, "CG", 1)
    F = Function(Q)
    F.interpolate(Expression(scalar_expression, degree=1))
    xdmf = XDMFFile("scalar_cg1_" + file_suffix + ".xdmf")
    xdmf.write_experimental(F, XDMFFile.Encoding_ASCII)

    Q = FunctionSpace(mesh, "CG", 2)
    F = Function(Q)
    F.interpolate(Expression(scalar_expression, degree=1))
    xdmf = XDMFFile("scalar_cg2_" + file_suffix + ".xdmf")
    xdmf.write_experimental(F, XDMFFile.Encoding_ASCII)

def output_vector_xdmf(mesh, vector_expression, file_suffix):
    Q = VectorFunctionSpace(mesh, "DG", 0)
    F = Function(Q)
    F.interpolate(Expression(vector_expression, degree=1))
    xdmf = XDMFFile("vector_dg0_" + file_suffix + ".xdmf")
    xdmf.write_experimental(F, XDMFFile.Encoding_ASCII)

    Q = VectorFunctionSpace(mesh, "DG", 1)
    F = Function(Q)
    F.interpolate(Expression(vector_expression, degree=1))
    xdmf = XDMFFile("vector_dg1_" + file_suffix + ".xdmf")
    xdmf.write_experimental(F, XDMFFile.Encoding_ASCII)

    Q = VectorFunctionSpace(mesh, "DG", 2)
    F = Function(Q)
    F.interpolate(Expression(vector_expression, degree=1))
    xdmf = XDMFFile("vector_dg2_" + file_suffix + ".xdmf")
    xdmf.write_experimental(F, XDMFFile.Encoding_ASCII)

    Q = VectorFunctionSpace(mesh, "CG", 1)
    F = Function(Q)
    F.interpolate(Expression(vector_expression, degree=1))
    xdmf = XDMFFile("vector_cg1_" + file_suffix + ".xdmf")
    xdmf.write_experimental(F, XDMFFile.Encoding_ASCII)

    Q = VectorFunctionSpace(mesh, "CG", 2)
    F = Function(Q)
    F.interpolate(Expression(vector_expression, degree=1))
    xdmf = XDMFFile("vector_cg2_" + file_suffix + ".xdmf")
    xdmf.write_experimental(F, XDMFFile.Encoding_ASCII)


# 2D
mesh = UnitSquareMesh(10, 10)
output_scalar_xdmf(mesh, "sin(2*pi*x[0])*cos(2*pi*x[1])", "2D")
output_vector_xdmf(mesh, ("sin(2*pi*x[0])","cos(2*pi*x[1])"), "2D")

# 2D manifold
mesh = UnitDiscMesh(mpi_comm_world(), 5, 1, 3)
output_scalar_xdmf(mesh, "sin(2*pi*x[0])*cos(2*pi*x[1])", "2D3")
output_vector_xdmf(mesh, ("sin(2*pi*x[0])","cos(2*pi*x[1])","1.0"), "2D3")

# 3D
mesh = UnitCubeMesh(5, 5, 5)
output_scalar_xdmf(mesh, "sin(2*pi*x[0])*cos(2*pi*x[1])", "3D")
output_vector_xdmf(mesh, ("sin(2*pi*x[0])","cos(2*pi*x[1])", "1.0"), "3D")
