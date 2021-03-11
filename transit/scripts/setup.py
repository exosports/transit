#!/usr/bin/env pyhton

from distutils.core import setup, Extension
import os.path as op
import glob
import numpy

transit_objs = [op.realpath(f) for f in glob.glob("./src/*.o")]
pu_objs = [op.realpath(f) for f in glob.glob("../pu/src/*.o")
      if not op.basename(f).startswith("messagep")]

transit_module = Extension('_transit_module',
      sources = ['python/transit_wrap.c'],
      include_dirs=[numpy.get_include()],
      extra_objects = transit_objs + pu_objs)

setup (name="transit_module",
       version= '0.1',
       author = 'Nate Lust',
       description = """A wrapper around the transit source code""",
       ext_modules = [transit_module],
       py_modules  = ["transit_module"],)
