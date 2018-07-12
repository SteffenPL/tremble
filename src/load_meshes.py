import fenics as fe

#
# class MuscleProblem:
#
#     def __init__(self):
#
#         self.mesh = None
#         self.dx = None
#         self.ds = None
#         self.boundaries = None
#
#         # boundary conditions
#         self.bc = None
#
#         self.Body_force = None
#
#
#     def add_boundary_condition(self):
#
#         pass



def load_2d_muscle_geo(filename='../geo/muscle_2d.xml', L0=1e-2):

    mesh = fe.Mesh(filename)
    coords = mesh.coordinates()
    coords *= L0
    mesh.bounding_box_tree().build(mesh)

    # Define Boundaries
    bottom = fe.CompiledSubDomain("near(x[1], side, 0.01) && on_boundary", side=-20.0 * L0)
    top = fe.CompiledSubDomain("near(x[1], side, 0.01) && on_boundary", side=20.0 * L0)

    # Initialize mesh function for boundary domains
    boundaries = fe.MeshFunction('size_t', mesh, 2)
    boundaries.set_all(0)
    bottom.mark(boundaries, 1)
    top.mark(boundaries, 2)

    # Define new measures associated with the interior domains and
    # exterior boundaries

    dx = fe.Measure('dx', domain=mesh)
    ds = fe.Measure('ds', domain=mesh, subdomain_data=boundaries)

    return mesh, dx, ds, {"top": top, "bottom": bottom}

def load_2d_muscle_bc(V_u, V_p, V_v, boundaries):

    bcs_u = []
    bcs_p = []
    bcs_v = []

    if V_u is not None:

        top_bc = fe.Expression(('0', '0'), element=V_u.ufl_element())
        bottom_bc = fe.Expression(('0', '0'), element=V_u.ufl_element())

        bcs_u.append(fe.DirichletBC(V_u, top_bc, boundaries['top']))
        bcs_u.append(fe.DirichletBC(V_u, bottom_bc, boundaries['bottom']))

    if V_v is not None:

        top_bc = fe.Expression(('0', '0'), element=V_v.ufl_element())
        bottom_bc = fe.Expression(('0', '0'), element=V_v.ufl_element())

        bcs_v.append(fe.DirichletBC(V_v, top_bc, boundaries['top']))
        bcs_v.append(fe.DirichletBC(V_v, bottom_bc, boundaries['bottom']))

    return (bcs_u, bcs_p, bcs_v)
