#!/usr/bin/env pyhton

from distutils.core import setup, Extension
import os.path as op
src = op.realpath(op.join(op.dirname(\
      op.realpath(__file__)),'../src'))
pu  = op.realpath(op.join(op.dirname(\
      op.realpath(__file__)),'../../pu/src'))


transit_module = Extension('_transit_module',sources=['src/transit_wrap.c'],
                            extra_objects=[op.join(src,'argum.o'),op.join(src,'cia.o'),
                                           op.join(src,'eclipse.o'),
                                           op.join(src,'extinction.o'),
                                           op.join(src,'geometry.o'),
                                           op.join(src,'idxrefraction.o'),
                                           op.join(src,'makesample.o'),
                                           op.join(src,'opacity.o'),
                                           op.join(src,'readatm.o'),
                                           op.join(src,'readlineinfo.o'),
                                           op.join(src,'slantpath.o'),
                                           op.join(src,'tau.o'),
                                           op.join(src,'transit.o'),
                                           op.join(src,'transitstd.o'),
                                           op.join(pu,'voigt.o'),
                                           op.join(pu,'procopt.o'),
                                           op.join(pu,'iomisc.o'),
                                           op.join(pu,'numerical.o'),
                                           op.join(pu,'xmalloc.o'),
                                           op.join(pu,'spline.o'),
                                           op.join(pu,'sampling.o')])

setup (name="transit_module",
       version= '0.1',
       author = 'Nate Lust',
       description = """A wrapper around the transit source code""",
       ext_modules = [transit_module],
       py_modules  = ["transit_module"],)
