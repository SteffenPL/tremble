import fenics as fe
from fenics import grad, inner, div, dot, tr, det
from math import pi
import numpy as np
import matplotlib.pyplot as plt
import time

from load_meshes import load_2d_muscle_geo, load_2d_muscle_bc

from plotting_tools import *



mesh, dx, ds, boundaries = load_2d_muscle_geo()
cell = mesh.ufl_cell()
dim = 2

# Notation
#  u: displacement field
#  p: lagrange multiplier, pressure
#  v: velocity field

element_u = fe.VectorElement('CG', cell, 2)
element_p = fe.FiniteElement('CG', cell, 1)
element_v = fe.VectorElement('CG', cell, 2)

# define mixed elements
mix_up  = fe.MixedElement([element_u, element_p           ])
mix_uv  = fe.MixedElement([element_u,            element_v])
mix_pv  = fe.MixedElement([           element_p, element_v])
mix_upv = fe.MixedElement([element_u, element_p, element_v])

# define function spaces
V_u = fe.FunctionSpace(mesh, element_u)
V_p = fe.FunctionSpace(mesh, element_p)
V_v = fe.FunctionSpace(mesh, element_v)
V_up = fe.FunctionSpace(mesh, mix_up)
V_uv = fe.FunctionSpace(mesh, mix_uv)
V_pv = fe.FunctionSpace(mesh, mix_pv)
V_upv = fe.FunctionSpace(mesh, mix_upv)

# Material constants
c_10 = fe.Constant(1e4)
c_01 = fe.Constant(2e3)

rho = fe.Constant(1.)

B = fe.Constant((0.,0.))


# bad class design, properties return copies!
class Solution(fe.Function):
    def __init__(self, u=None, p=None, v=None, t=None):
        super().__init__(V_upv)
        self.reset()

        if u is not None:
            self.u = u
        if p is not None:
            self.p = p

        if v is not None:
            self.v = v

        if t is not None:
            self.t = t


    def reset(self):
        self.__u = None
        self.__p = None
        self.__v = None
        self.__up = None
        self.__pv = None
        self.__uv = None
        self.t = 0
        self.__upv = self

    @property
    def u(self):
        if self.__u is None and self.__upv is not None:
            self.__u = fe.Function(V_u)
        fe.assign(self.__u, self.__upv.sub(0))
        return self.__u

    @property
    def p(self):
        if self.__p is None and self.__upv is not None:
            self.__p = fe.Function(V_p)
        fe.assign(self.__p, self.__upv.sub(1))
        return self.__p

    @property
    def v(self):
        if self.__v is None and self.__upv is not None:
            self.__v = fe.Function(V_v)
        fe.assign(self.__v, self.__upv.sub(2))
        return self.__v

    @property
    def up(self):
        if self.__up is None and self.__upv is not None:
            self.__up = fe.Function(V_up)
        fe.assign(self.__up.sub(0), self.__upv.sub(0))
        fe.assign(self.__up.sub(1), self.__upv.sub(1))
        return self.__up

    @property
    def uv(self):
        if self.__uv is None and self.__upv is not None:
            self.__uv = fe.Function(V_uv)
        fe.assign(self.__uv.sub(0), self.__upv.sub(0))
        fe.assign(self.__uv.sub(1), self.__upv.sub(2))
        return self.__uv


    @property
    def pv(self):
        if self.__pv is None and self.__upv is not None:
            self.__pv = fe.Function(V_pv)

        fe.assign(self.__pv.sub(0), self.upv.sub(1))
        fe.assign(self.__pv.sub(1), self.upv.sub(2))
        return self.__pv

    @property
    def upv(self):
        return self.__upv

    @upv.setter
    def upv(self, val):
        self.reset()
        fe.assign(self.__upv, val)

    @u.setter
    def u(self, val):
        self.__u = None
        self.__up = None
        self.__uv = None

        fe.assign(self.__upv.sub(0), val)

    @p.setter
    def p(self, val):
        self.__p = None
        self.__up = None
        self.__pv = None

        fe.assign(self.__upv.sub(1), val)

    @v.setter
    def v(self, val):
        self.__v = None
        self.__uv = None
        self.__pv = None

        fe.assign(self.__upv.sub(2), val)

    @up.setter
    def up(self, val):
        self.u = val.sub(0)
        self.p = val.sub(1)

    @uv.setter
    def uv(self, val):
        self.u = val.sub(0)
        self.v = val.sub(1)

    @pv.setter
    def pv(self, val):
        self.p = val.sub(0)
        self.v = val.sub(1)




def deformation_grad(u):
    I = fe.Identity(dim)
    return fe.variable(I + grad(u))


def invariants(F):
    C = F.T * F

    # Invariants
    I_1 = tr(C) + (3-dim)
    I_2 = 0.5*((tr(C)+(3-dim))**2 - (tr(C*C)+(3-dim)))
    J = det(F)  # this is a third invariant 'I_3'

    return I_1, I_2, J

def isochronic_deformation_grad(F, J):
    return J**(-1./3.)*F


def material_mooney_rivlin(I_1, I_2, c_10, c_01):
    # Mooney-Rivlin material law
    return c_10*(I_1 - 3) + c_01*(I_2 - 3)


def kinetic_enegery(v, rho):
    return 0.5 * rho * inner(v, v)


def incompr_constr(J):
    return (J-1)


def incompr_relaxation(p, kappa):
    return 0.5*p**2/kappa


def first_piola_stress(L, F):
    return -fe.diff(L, F)


def incompr_stress(g, F):
    return fe.diff(g, F)


def weak_div_term(P, eta):
    return -inner(P, grad(eta))*dx


def const_eq(L_mod, p):
    return fe.diff(L_mod, p)


'''
    We now collect boundary conditions
'''



'''
    Now we have define the Lagrange function. Then we are ready to generate
    the weak for of the equations. For convenience, we use the notation of elasticity.
'''

def teast_G():
    assert dim==2, 'this test only works in two dimensions.'


def solve_steady_state_heiser_weissinger(kappa):

    w = fe.Function(V_up)
    (u, p) = fe.split(w)
    p = fe.variable(p)
    (eta, q) = fe.TestFunctions(V_up)
    dw = fe.TrialFunction(V_up)

    kappa = fe.Constant(kappa)

    bcs_u, bcs_p, bcs_v = load_2d_muscle_bc(V_up.sub(0), V_up.sub(1), None, boundaries)

    F = deformation_grad(u)
    I_1, I_2, J = invariants(F)
    F_iso = isochronic_deformation_grad(F, J)
    I_1_iso, I_2_iso  = invariants(F)[0:2]

    W = material_mooney_rivlin(I_1_iso, I_2_iso, c_10, c_01) #+ incompr_relaxation(p, kappa)
    g = incompr_constr(J)

    # Lagrange function (without constraint)
    L = -W

    # Modified Lagrange function (with constraints)
    L_mod = L - p*g
    P = first_piola_stress(L, F)
    G = incompr_stress(g, F)  # = J*fe.inv(F.T)

    Lp = const_eq(L_mod, p)

    a_static = weak_div_term(P + p*G, eta) + inner(B, eta)*dx + inner(Lp, q)*dx

    J_static = fe.derivative(a_static, w, dw)
    ffc_options = {"optimize": True}
    problem = fe.NonlinearVariationalProblem(a_static, w, bcs_u + bcs_p, J=J_static, form_compiler_parameters=ffc_options)
    solver = fe.NonlinearVariationalSolver(problem)

    solver.solve()

    return w



def explicit_relax_dyn(w0, kappa=1e5, dt=1.e-5, t_end=1.e-4, show_plots=False):

    (u0, p0, v0) = fe.split(w0)

    bcs_u, bcs_p, bcs_v = load_2d_muscle_bc(V_upv.sub(0), V_upv.sub(1), V_upv.sub(2), boundaries)

    kappa = fe.Constant(kappa)


    F = deformation_grad(u0)
    I_1, I_2, J = invariants(F)
    F_iso = isochronic_deformation_grad(F, J)
    #I_1_iso, I_2_iso  = invariants(F_iso)[0:2]

    W = material_mooney_rivlin(I_1, I_2, c_10, c_01) + incompr_relaxation(p0, kappa)
    g = incompr_constr(J)

    # Lagrange function (without constraint)

    L = -W
    P = first_piola_stress(L, F)
    G = incompr_stress(g, F)

    (u1, p1, v1) = fe.TrialFunctions(V_upv)
    (eta, q, xi) = fe.TestFunctions(V_upv)
    a_dyn_u = inner(u1-u0, eta)*dx - dt*inner(v1, eta)*dx
    a_dyn_p = (p1-p0)*q*dx - dt*kappa*div(v1)*J*q*dx
    #a_dyn_v = rho*inner(v1-v0, xi)*dx + dt*(inner(P + p0*G, grad(xi))*dx - inner(B, xi)*dx)
    a_dyn_v = rho*inner(v1-v0, xi)*dx + dt*(inner(P, grad(xi))*dx + inner(p0*G,grad(xi))*dx - inner(B, xi)*dx)



    a = fe.lhs(a_dyn_u + a_dyn_p + a_dyn_v)
    l = fe.rhs(a_dyn_u + a_dyn_p + a_dyn_v)

    w1 = fe.Function(V_upv)
    w2 = fe.Function(V_upv)

    sol = []

    vol = fe.assemble(1.*dx)


    t = 0
    while t < t_end:
        print("progress: %f" % (100.*t/t_end))

        A = fe.assemble(a)
        L = fe.assemble(l)

        for bc in bcs_u + bcs_p + bcs_v:
            bc.apply(A, L)

        fe.solve(A, w1.vector(), L)

        if fe.norm(w1.vector()) > 1e7:
            print('ERROR: norm explosion')
            break

        # update initial values for next step
        w0.assign(w1)
        t += dt

        if show_plots:
            # plot result
            fe.plot(w0.sub(0), mode='displacement')
            plt.show()

        # save solutions
        sol.append(Solution(t=t))
        sol[-1].upv.assign(w0)

    return sol, W, kappa


def impl_dyn(w0, dt=1.e-5, t_end=1.e-4, show_plots=False):

    (u0, p0, v0) = fe.split(w0)

    bcs_u, bcs_p, bcs_v = load_2d_muscle_bc(V_upv.sub(0), V_upv.sub(1), V_upv.sub(2), boundaries)


    # Lagrange function (without constraint)

    (u1, p1, v1) = fe.TrialFunctions(V_upv)
    (eta, q, xi) = fe.TestFunctions(V_upv)

    F = deformation_grad(u1)
    I_1, I_2, J = invariants(F)
    F_iso = isochronic_deformation_grad(F, J)
    #I_1_iso, I_2_iso  = invariants(F_iso)[0:2]

    W = material_mooney_rivlin(I_1, I_2, c_10, c_01)
    g = incompr_constr(J)
    L = -W
    P = first_piola_stress(L, F)
    G = incompr_stress(g, F)


    a_dyn_u = inner(u1-u0, eta)*dx - dt*inner(v1, eta)*dx
    a_dyn_p = inner(g,q)*dx
    a_dyn_v = rho*inner(v1-v0, xi)*dx + dt*(inner(P, grad(xi))*dx + inner(p1*G,grad(xi))*dx - inner(B, xi)*dx)

    a = a_dyn_u + a_dyn_p + a_dyn_v

    w1 = fe.Function(V_upv)

    sol = []


    t = 0
    while t < t_end:
        print("progress: %f" % (100.*t/t_end))

        fe.solve(a==0, w1, bcs_u + bcs_p + bcs_v)

        if fe.norm(w1.vector()) > 1e7:
            print('ERROR: norm explosion')
            break

        # update initial values for next step
        w0.assign(w1)
        t += dt

        if show_plots:
            # plot result
            fe.plot(w0.sub(0), mode='displacement')
            plt.show()

        # save solutions
        sol.append(Solution(t=t))
        sol[-1].upv.assign(w0)

    return sol, W, kappa



def half_exp_dyn(w0, dt=1.e-5, t_end=1.e-4, show_plots=False):

    u0 = w0.u
    p0 = w0.p
    v0 = w0.v

    bcs_u, bcs_p, bcs_v = load_2d_muscle_bc(V_u, V_pv.sub(0), V_pv.sub(1), boundaries)

    F = deformation_grad(u0)
    I_1, I_2, J = invariants(F)
    F_iso = isochronic_deformation_grad(F, J)
    #I_1_iso, I_2_iso  = invariants(F_iso)[0:2]
    W = material_mooney_rivlin(I_1, I_2, c_10, c_01)
    g = incompr_constr(J)
    # Lagrange function (without constraint)
    L = -W
    P = first_piola_stress(L, F)
    G = incompr_stress(g, F)

    # a_dyn_u = inner(u1 - u0, eta) * dx - dt * inner(v1, eta) * dx


    u1 = fe.TrialFunction(V_u)
    eta = fe.TestFunction(V_u)

    #u11 = fe.Function(V_u)
    #F1 = deformation_grad(u11)
    #g1 = incompr_constr(fe.det(F1))
    #G1 = incompr_stress(g1, F1)


    (p1, v1) = fe.TrialFunctions(V_pv)
    (q, xi) = fe.TestFunctions(V_pv)

    a_dyn_u = inner(u1-u0,eta)*dx - dt*inner(v0,eta)*dx

    a_dyn_p = fe.tr(G*grad(v1))*q*dx
    #a_dyn_v = rho*inner(v1-v0, xi)*dx + dt*(inner(P + p0*G, grad(xi))*dx - inner(B, xi)*dx)
    a_dyn_v = rho*inner(v1-v0, xi)*dx + dt*(inner(P, grad(xi))*dx + inner(p1*G,grad(xi))*dx - inner(B, xi)*dx)


    a_u = fe.lhs(a_dyn_u)
    l_u = fe.rhs(a_dyn_u)

    a_pv = fe.lhs(a_dyn_p + a_dyn_v)
    l_pv = fe.rhs(a_dyn_p + a_dyn_v)


    u1 = fe.Function(V_u)
    pv1 = fe.Function(V_pv)

    sol = []

    vol = fe.assemble(1.*dx)

    A_u = fe.assemble(a_u)
    A_pv = fe.assemble(a_pv)

    for bc in bcs_u:
        bc.apply(A_u)

    t = 0
    while t < t_end:
        print("progress: %f" % (100.*t/t_end))

        # update displacement u
        L_u = fe.assemble(l_u)
        fe.solve(A_u, u1.vector(), L_u)
        u0.assign(u1)


        L_pv = fe.assemble(l_pv)
        for bc in bcs_p + bcs_v:
            bc.apply(A_pv, L_pv)

        fe.solve(A_pv, pv1.vector(), L_pv)

        if fe.norm(pv1.vector()) > 1e8:
            print('ERROR: norm explosion')
            break

        # update initial values for next step
        w0.u = u1
        w0.pv = pv1
        p0.assign(w0.p)
        v0.assign(w0.v)


        t += dt

        if show_plots:
            # plot result
            fe.plot(w0.sub(0), mode='displacement')
            plt.show()

        # save solutions
        sol.append(Solution(t=t))
        sol[-1].upv.assign(w0)

    return sol, W



B.assign(fe.Constant((10000, 0)))

w = solve_steady_state_heiser_weissinger(kappa=1e7)
fe.plot(w.sub(0), mode='displacement')
plt.show()


B.assign(fe.Constant((0,0)))

w0 = Solution()
w0.up = w
w0.t = 0

dt = 2.e-5
sol, W = half_exp_dyn(w0, dt=dt, t_end=500*dt, show_plots=False)
# dif = explicit_relax_dyn(w0, kappa=1e4, dt=dt, t_end=1000*dt, show_plots=False)

(u0, p0, v0) = fe.split(w0)
T = 0.5*rho*inner(v0,v0)

plot_form(sol, w0, W*dx)
plot_form(sol, w0, T*dx)
plot_form(sol, w0, (T+W)*dx)
plot_form(sol, w0, (det(deformation_grad(u0)))*dx)



save_to_vtk(sol)

"""
# def solve_dynamics():

u0 =  w.split(deepcopy=True)[0]
# p = w.sub(1, deepcopy=True)
v0 = fe.Function(V_v)


(p1, v1) = fe.TrialFunctions(V_pv)
(q, xi) = fe.TestFunctions(V_pv)

bcs_u, bcs_p, bcs_v = load_2d_muscle_bc(V_u, V_pv.sub(0), V_pv.sub(1), boundaries)



dt = 0.001
t_end = 1


w1 = fe.Function(V_pv)
u1 = fe.Function(V_u)

B.assign(fe.Constant((0,0)))

bcs_u, bcs_p, bcs_v = load_2d_muscle_bc(V_up.sub(0), V_up.sub(1), None, boundaries)

sol = []
t = 0
while t < 10*dt:

    F0 = deformation_grad(u0)
    I_10, I_20, J0 = invariants(F0)
    W0 = material_mooney_rivlin(I_10, I_20, c_10, c_01)
    g0 = incompr_constr(J0)

    P0 = first_piola_stress(W0, F0)
    G0 = incompr_stress(g0, F0)

    # update position
    u1.assign(u0)
    u1.vector().axpy(dt, v0.vector())

    # TODO: either u, v are discretizised in the same way or interpolate needed!

    F1 = deformation_grad(u1)
    J1 = fe.det(F1)
    g1 = incompr_constr(J1)

    G1 = incompr_stress(g1, F1)



    # update velocity and lagrange multiplier
    # a_v = inner(v1, xi)*dx + inner(dt*p1*G0, grad(xi))*dx
    # a_p = p1*q*dx  # inner(G1, grad(v1))*q*dx

    # l_v = inner(v0, xi)*dx + dt*weak_div_term(P0, xi) + dt*inner(B, xi)*dx

    # fe.solve( a_v+a_p == l_v, w1, bcs_p + bcs_v)

    # U1 = fe.Function(V_u)
    # Y1 = fe.Function(V_pv)
    # Y1.assign(u0)
    # Y1.assign(w1)
    # sol.append((U1, Y1, t))

    # assign new values as next initial values
    u0.assign(u1)
    v0.assign(w1.sub(1, deepcopy=True))

    fe.plot(u1, mode='displacement')
    plt.show()

    t += dt

def save_to_vtk(sol):
    fu = fe.File('/tmp/u.pvd')
    fp = fe.File('/tmp/p.pvd')
    fv = fe.File('/tmp/v.pvd')

    for U, Y, t in sol:
        fu << U, t
        fp << Y.sub(0), t
        fv << Y.sub(1), t

save_to_vtk(sol)


def plt_show(i):
    fe.plot(sol[i][0], mode='displacement')
    plt.show()


plt_show(0)

"""
