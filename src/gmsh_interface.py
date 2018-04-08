import os
import re
import tempfile
import subprocess
import shlex


class GmshInterface:

    def __init__(self, fn: str, dim: int = 3, dolfin_convert_exec: str = None):
        assert ".geo" in fn, "Filetype .geo expected."

        self.dim = dim
        self.fn = fn
        self.parameters = dict()
        self.pattern =  re.compile(r"\s*(\S*)\s* = \s*DefineNumber\[\s*(\S*)\s*," , re.VERBOSE )

        self.find_parameters()

        if dolfin_convert_exec is None:
            self.dolfin_convert_exec = "/opt/miniconda3/envs/fenics/bin/dolfin-convert"


    def find_parameters(self):

        f = open(self.fn, "r")

        for line in f:
            m = self.pattern.match(line)
            if m is not None:
                key = m.groups()[0]
                value = m.groups()[1]

                self.set_parameter(key, value)

    def set_parameter(self, key: str, value: float):
        self.parameters[key] = value

    def get_parameters(self, key: str) -> float:
        return self.parameters[key]

    def generate_msh(self, fn_output: str = None, lc: float = None):

        if fn_output is None:
            fn_output = self.fn.replace(".geo",".msh")

        if ".msh" not in fn_output:
            fn_output = fn_output + ".msh"

        if lc is not None:
            self.set_parameter("lc",lc)

        # write modified .geo file

        fin = open(self.fn,"r")
        f_tmp = tempfile.NamedTemporaryFile("w",delete=False)

        for line in fin:
            m = self.pattern.match(line)
            if m is not None:
                key = m.groups()[0]

                assert key in self.parameters, "Parameter not detected on first parse."

                line = key + "=" + str(self.parameters[key]) + ";\n"
            f_tmp.write(line)

        f_tmp.close()
        subprocess.call(["gmsh","-" + str(self.dim),f_tmp.name,"-o",fn_output])

    def generate_xml(self, fn_output: str = None,  lc: float = None):
        if ".xml" not in fn_output:
            fn_output = fn_output + ".xml"

        fn_msh = fn_output.replace(".xml",".msh")
        self.generate_msh(fn_msh)

        subprocess.call([self.dolfin_convert_exec,fn_msh,fn_output])






