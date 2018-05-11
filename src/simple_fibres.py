import fenics as fe
# We load a few fenics objects into the global namespace
from fenics import div, grad, curl, inner, dot, inv, tr
import numpy as np

#import matplotlib.pyplot as plt

# some fenics settings
fe.parameters['form_compiler']['cpp_optimize'] = True
fe.parameters['form_compiler']['optimize'] = True
#fe.parameters["form_compiler"]["quadrature_degree"] = 5


# Create mesh and define function space
n = 20
mesh = fe.UnitCubeMesh(n, n, n)

# Init function spaces
element_3 = fe.VectorElement("P", mesh.ufl_cell(), 1)
element = fe.FiniteElement("P", mesh.ufl_cell(), 2)

# Mixed function space
TH = element_3 * element
V = fe.FunctionSpace(mesh, TH)

# Define Boundaries
left =  fe.CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = fe.CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)

# Define Dirichlet boundary (x = 0 or x = 1)
u_left = fe.Expression(("0.0", "0.0", "0.0"), element=element_3)
u_right = fe.Expression(("0.0", "0.0", "0.0"), element=element_3)
p_left = fe.Constant(0.)

# Define acting force
b = fe.Constant((0.0, 0.0, 0.0))       # Body force per unit volume
t_bar  = fe.Constant((0.0,  0.0, 0.0))   # Traction force on the boundary


# Define test and trial functions
w = fe.Function(V)  # most recently computed solution
(u,p) = fe.split(w)
(v,q) = fe.TestFunctions(V)
dw = fe.TrialFunction(V)


# Kinematics
d = u.geometric_dimension()
I = fe.Identity(d)             # Identity tensor

F = fe.variable(I + grad(u))  # Deformation gradient
C = fe.variable(F.T*F)  # Right Cauchy-Green tensor
J = fe.det(C)

DE = lambda v: 0.5*(F.T*grad(v) + grad(v).T*F)

a_0 = fe.as_vector([[1.0], [0.], [0.]])

# Invariants
I_1 = tr(C)
I_2 = 0.5*(tr(C)**2 - tr(C*C))

I_4 = dot(a_0.T, C*a_0)

# Continuation parameter
lmbda = 0.01

# parameters taken from O. Roehrle and Heidlauf

# Mooney-Rivlin parameters (in kPa)

b_1 = 2.756e-5
d_1 = 43.373

p_max = 73.0
f_s = 50.  # stimulation frequency (in Hz)

# Not found
rho = fe.Constant(1.1e3)
c_10 = 6.352e-10
c_01 = 3.627
kappa = fe.Constant(0.01)

# currently simplifications
p_act = 0.
lmb_f = 2.

W_iso = c_10*(I_1 - 3) + c_01*(I_2 - 3) + 0.5*kappa*(J-1)**2

dx = fe.dx
ds = fe.ds

P = fe.diff(W_iso,F)
S = 2 * fe.diff(W_iso,C)
F_static = inner(S, DE(v))*dx - rho*inner(b, v)*dx - inner(t_bar, v)*ds + (p/kappa + J - 1)*q*dx
#F_static = inner(P, grad(v))*dx - rho*inner(b, v)*dx - inner(t_bar, v)*ds + (p/kappa + J - 1)*q*dx

file_u = fe.File("../results/displacement.pvd")
file_p = fe.File("../results/pressure.pvd")

bcul = fe.DirichletBC(V.sub(0), u_left, left)
bcur = fe.DirichletBC(V.sub(0), u_right, right)
#bcpl = fe.DirichletBC(V.sub(1), p_left, left)
bcs = [bcul, bcur]  # bcpl]

J_static = fe.derivative(F_static,w,dw)

#fe.solve(F_static==0, w, bcs)

ffc_options = {"optimize": True}
problem = fe.NonlinearVariationalProblem(F_static, w, bcs, J=J_static, form_compiler_parameters=ffc_options)
solver = fe.NonlinearVariationalSolver(problem)

# set parameters
prm = solver.parameters
if True:
    iterative_solver = True
    #prm['newton_solver']['absolute_tolerance'] = 1E-3
    #prm['newton_solver']['relative_tolerance'] = 1E-2
    #prm['newton_solver']['maximum_iterations'] = 20
    #prm['newton_solver']['relaxation_parameter'] = 1.0
    prm['newton_solver']['linear_solver'] = 'bicgstab'
    #if iterative_solver:
    #    prm['newton_solver']['krylov_solver']['absolute_tolerance'] = 1E-4
    #    prm['newton_solver']['krylov_solver']['relative_tolerance'] = 1E-5
    #    prm['newton_solver']['krylov_solver']['maximum_iterations'] = 1000

fe.set_log_level(fe.PROGRESS)

solver.solve()

(u,p) = w.split()

file_u << u
file_p << p