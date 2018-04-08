from gmsh_interface import GmshInterface

g = GmshInterface("../geo/simple_muscle_3d.geo")
g.set_parameter("lc",0.1)
g.generate_xml("../geo/test.xml")




from fenics import Mesh, plot
import matplotlib.pyplot as plt

mesh = Mesh("../geo/test.xml")
plot(mesh)
plt.show()