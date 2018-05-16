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
n = 10
mesh = fe.UnitCubeMesh(n, n, n)

# Init function spaces
element_3 = fe.VectorElement("P", mesh.ufl_cell(), 1)
element = fe.FiniteElement("DP", mesh.ufl_cell(), 0)

# Mixed function space
TH = element_3 * element
V = fe.FunctionSpace(mesh, TH)

# Define Boundaries
left =  fe.CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = fe.CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)

# Define acting force
b = fe.Expression(("scale*gravity", "0.", "0."), scale=1., gravity=-0, element=element_3)       # Body force per unit volume
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
J = fe.det(F)

DE = lambda v: 0.5*(F.T*grad(v) + grad(v).T*F)

a_0 = fe.as_vector([[1.0], [0.], [0.]])
A_00 = dot( a_0, a_0.T )

# Invariants
I_1 = tr(C)
I_2 = 0.5*(tr(C)**2 - tr(C*C))

I_4 = dot(a_0.T, C*a_0)

# parameters taken from O. Roehrle and Heidlauf

# Mooney-Rivlin parameters (in Pa)
#c_10 = fe.Constant(6.352e4)
#c_01 = fe.Constant(3.627e3)
c_10 = fe.Constant(1e4)
c_01 = fe.Constant(2e2)
b_1 = 2.756e-5
d_1 = 43.373
p_max = 73.0
f_s = 50.  # stimulation frequency (in Hz)

# Not found
rho = fe.Constant(1.1e3)
kappa = fe.Expression("kappa", kappa=1e8, degree=0)

# currently simplifications
p_act = fe.Expression("p_act*scale", p_act=1., scale=1., degree=0)
lmb_f = fe.Expression("lmd_f", lmd_f=2.0, scale=1., degree=0)

I1_ = J**(-2./3.)*I_1
I2_ = J**(-4./3.)*I_2

W_iso = c_10*(I1_ - 3) + c_01*(I2_ - 3) - p*(J-1) - p**2/2*kappa


dx = fe.dx
ds = fe.ds

P_iso = fe.diff(W_iso, F)
S_iso = F*P_iso

S_passive = lmb_f**(-1)*p_act*A_00


S = S_iso + S_passive

F_static = inner(S, DE(v))*dx - rho*inner(b, v)*dx - inner(t_bar, v)*ds + (p/kappa + (J-1))*q*dx

file_u = fe.File("../results/displacement.pvd")
file_p = fe.File("../results/pressure.pvd")

# Define Dirichlet boundary (x = 0 or x = 1)
u_left = fe.Expression(("0.0", "0.0", "0.0"), element=element_3)
u_right = fe.Expression(("1.*scale", "0.0", "0.0"), scale=0.2, element=element_3)
p_left = fe.Constant(0.)

bcul = fe.DirichletBC(V.sub(0), u_left, left)
bcur = fe.DirichletBC(V.sub(0), u_right, right)
bcs = [bcul, ] #bcur]

J_static = fe.derivative(F_static,w,dw)


ffc_options = {"optimize": True}
problem = fe.NonlinearVariationalProblem(F_static, w, bcs, J=J_static, form_compiler_parameters=ffc_options)
solver = fe.NonlinearVariationalSolver(problem)

# set parameters
prm = solver.parameters
if True:
    iterative_solver = True
    #prm['newton_solver']['absolute_tolerance'] = 1e-10
    #prm['newton_solver']['relative_tolerance'] = 1e-9
    #prm['newton_solver']['maximum_iterations'] = 20
    #prm['newton_solver']['relaxation_parameter'] = 1.0
    prm['newton_solver']['linear_solver'] = 'bicgstab'
    #if iterative_solver:
    #    prm['newton_solver']['krylov_solver']['absolute_tolerance'] = 1E-4
    #    prm['newton_solver']['krylov_solver']['relative_tolerance'] = 1E-5
    #    prm['newton_solver']['krylov_solver']['maximum_iterations'] = 1000

fe.set_log_level(fe.PROGRESS)


prog = 0.

dprog = 0.2

w_last = fe.Function(V)

while prog < 1:

    success = False


    # update parameters
    while success is False:
        kappa.kappa = 1e6 + 0*(prog+dprog)*1e3
        u_right.scale = (prog+dprog)
        #b.scale = (prog+dprog)
        p_act.scale = 100*(prog+dprog)

        print('dprop',dprog, '  prog ', prog)

        try:
            solver.solve()
            success = True
            w_last.assign(w)
            prog += dprog

            (u, p) = w.split()

            file_u << u, prog
            file_p << p, prog

            dprog *= 2
            dprog = min(dprog, 1.-prog)

        except RuntimeError:
            print("Solve failed, try smaller stepsize")
            dprog /= 2
            w.assign(w_last)



file_J = fe.File("../results/det.pvd")
file_J << fe.project(J, fe.FunctionSpace(mesh, "P", 1)), prog


file_E = fe.File("../results/strain_energy.pvd")
file_E << fe.project(W_iso, fe.FunctionSpace(mesh, "P", 1))


file_S = fe.File("../results/stress.pvd")
file_S << fe.project(S, fe.TensorFunctionSpace(mesh, "DP", 0))
