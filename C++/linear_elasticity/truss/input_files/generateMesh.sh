#!/bin/bash
rm 3D-truss.msh
gmsh 3D-truss.geo -o 3D-truss.msh -order 1 -3
