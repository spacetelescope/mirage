#! /usr/bin/env python

from mirage import imaging_simulator
t1 = imaging_simulator.ImgSim()
t1.paramfile = 'yaml_files/jw00793001001_01101_00001_nis.yaml'
#t1.paramfile = 'yaml_files/jw00793001001_01102_00001_nis.yaml'
t1.create()
