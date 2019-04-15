#!/bin/sh

./generateCSV.py > kinamatics_fem.csv
cp kinamatics_fem.csv /root/workspace/partons-example/data/examples/cff/kinematics_dvcs_cff.csv
cp examples.cpp /root/workspace/partons-example/src/.
cp main.cpp /root/workspace/partons-example/src/.
