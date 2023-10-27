#!/usr/bin/env python

from distutils.core import setup

setup(
        name         = "pyband",
        version      = "1.0",
        description  = "band plot using python matplotlib",
        author       = "Qijing Zheng",
        author_email = "zqj.kaka@gmail.com",
        url          = 'https://github.com/QijingZheng/VaspBandUnfolding',
        py_modules   = [],
        scripts      = [
            "aseconv.py",
            "energy_unit_conv.py",
            "kut.py",
            "molAdd.py",
            "molDel.py",
            "nebplt.py",
            "nose_mass.py",
            "plot_workfunc.py",
            "velInit.py",
            "xcell.py",
            "xtraj.py",
            "npband",
            "npdos",
            "pyband",
            "pydos",
            "pygap",
            ]
        )
