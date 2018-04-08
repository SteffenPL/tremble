#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 12:31:18 2018

@author: plunder
"""
import subprocess

import fenics as fe
from dolfin import Function
from fenics import near, between, Function, TestFunction, TrialFunction, grad, div, det, tr, dot

from gmsh_interface import GmshInterface


def convert_msh_to_xml(fn_msh, fn_xml):
    subprocess.call("dolfin-convert ../geo/%s ../geo/%s" % (fn_msh, fn_xml))


class Hyperplane(fe.SubDomain):
    def __init__(self,z_value):
        self.z_value = z_value

    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], self.z_value)


class Stripe(fe.SubDomain):
    def __init__(self,z_interval):
        self.z_interval = z_interval

    def inside(self, x, on_boundary):
        return on_boundary and between(x[2], self.z_interval[0], self.z_interval[1])


class MuscleState:

    def __init__(self):
        self.mesh = None

        # self.domains = None
        # self.boundaries = None

        self.dx = None
        self.ds = None

        # current deformation
        self.u = None

    def load_mesh(self, fn: str = None):
        if fn is not None:
            self.mesh = fe.Mesh(fn)
        else:
            self.mesh = fe.UnitCubeMesh(10,12,14)



    def compute_static_deformation(self):

        assert self.mesh is not None

        # now we define subdomains on the mesh

        bottom = fe.CompiledSubDomain('near(x[2], 0) && on_boundary')
        top = fe.CompiledSubDomain('near(x[2], 1) && on_boundary')
        # middle = fe.CompiledSubDomain('x[2] > 0.3 && x[2] < 0.7')

        # Initialize mesh function for interior domains
        self.domains = fe.MeshFunction('size_t', self.mesh, 3)
        self.domains.set_all(0)
        # middle.mark(self.domains, 1)

        # Initialize mesh function for boundary domains
        self.boundaries = fe.MeshFunction('size_t', self.mesh, 2)
        self.boundaries.set_all(0)
        bottom.mark(self.boundaries, 1)
        top.mark(self.boundaries, 2)

        # Define new measures associated with the interior domains and
        # exterior boundaries

        self.dx = fe.Measure('dx', domain=self.mesh, subdomain_data=self.domains)
        self.ds = fe.Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

        # define function spaces
        V = fe.VectorFunctionSpace(self.mesh, "Lagrange", 1)
        # now we define subdomains on the mesh

        bottom = fe.CompiledSubDomain('near(x[2], 0) && on_boundary')
        top = fe.CompiledSubDomain('near(x[2], 1) && on_boundary')
        # middle = fe.CompiledSubDomain('x[2] > 0.3 && x[2] < 0.7')

        d = self.mesh.geometry().dim()

        # Initialize mesh function for interior domains
        self.domains = fe.MeshFunction('size_t', self.mesh, d)
        self.domains.set_all(0)
        # middle.mark(self.domains, 1)

        # Initialize mesh function for boundary domains
        self.boundaries = fe.MeshFunction('size_t', self.mesh, d - 1)
        self.boundaries.set_all(0)
        bottom.mark(self.boundaries, 1)
        top.mark(self.boundaries, 2)

        # Define new measures associated with the interior domains and
        # exterior boundaries

        self.dx = fe.Measure('dx', domain=self.mesh, subdomain_data=self.domains)
        self.ds = fe.Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

        c_zero = fe.Constant((0,0,0))

        # define boundary conditions
        bc_bottom = fe.DirichletBC(V, c_zero, bottom)
        bc_top    = fe.DirichletBC(V, c_zero, top)

        bcs = [bc_bottom]  # , bc_top]

        # define functions
        du = TrialFunction(V)
        v = TestFunction(V)
        u = Function(V)
        B = fe.Constant((0., 2.0, 0.))
        T = fe.Constant((0.0, 0.0, 0.0))

        d = u.geometric_dimension()
        I = fe.Identity(d)
        F = I + grad(u)
        C = F.T*F

        I_1 = tr(C)
        J = det(F)

        E, mu = 10., 0.3
        mu, lmbda = fe.Constant(E/(2*(1+mu))), fe.Constant(E*mu/((1+mu)*(1-2*mu)))

        # stored energy (comp. neo-hookean model)
        psi = (mu/2.)*(I_1 - 3) - mu*fe.ln(J) + (lmbda/2.)*(fe.ln(J))**2

        dx = self.dx
        ds = self.ds

        Pi = psi*fe.dx - dot(B, u)*fe.dx - dot(T, u)*fe.ds

        F = fe.derivative(Pi, u, v)
        J = fe.derivative(F, u, du)

        fe.solve(F == 0, u, bcs, J=J)

        # save results
        self.u = u
        
        # write to disk
        output = fe.File("/tmp/static.pvd")
        output << u


g = GmshInterface("../geo/simple_muscle_3d.geo")
g.set_parameter("lc",0.08)
#for i in [1,2,3]:
#    g.set_parameter("T_"+str(i),0.4)
#    g.set_parameter("W_"+str(i),0.4)
g.generate_xml("../geo/test.xml")

import matplotlib.pyplot as plt


m = MuscleState()
m.load_mesh("../geo/test.xml")
fe.plot(m.mesh)
plt.show()
m.compute_static_deformation()
