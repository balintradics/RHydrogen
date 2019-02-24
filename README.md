# RHydrogen
RHydrogen is a lightweight open-source c++ program which solves the Dirac equation for Hydrogen-like atoms. 

## General:

This program solves the radial Dirac-equation with a single electron and a
Coulomb potential using the mapped Fourier-grid method. 
It is based on this paper:

J.Phys.A: Math.Gen. 38 (2005) 3157-3171

The diagonalization of the Dirac-Coulomb Hamiltonian is done ROOT (to be changed to GSL).


## Requirements:

ROOT (if possible >= 5.18), http://root.cern.ch

## Build:

(here I assume that the ROOT environment variables are set)

```
make
```

## Run:
```
./rhydrogen <parameters>
```
for help, type:
```
./rhydrogen --help
```
After finished the code creates a ROOT file (eigenvectors.root) in which 2x10 TGraphs are saved. 
These are the real part of the numerical bigger spinor eigenvector values and their squared as a 
function of the grid position (in hbar/mc), for the first 10 eigenvectors.

## Example session:

- Number of grid points is N = 50, 
- scale parameter is s = 4000 (consult the paper), 
- orbital angular momentum quantum number l = 0 (that is s orbit), 
- the 1/2 spin angular momentum quantum number subtracts off from l (j = l + (-1/2) = 0.5),  
- the element is Hydrogen (Z = 1).
```
./rhydrogen -N 50 -s 4000 -l 0 -j 0.5 -Z 1
Scale parameter: 4000
Grid points: 50
Orbital angular mom. quantum number: 0
Total angular mom. quantum number: 0.5
Z nucl. charge: 1
State: n s^{1/2}
Electron mass: 5.109989E+05 eV

=================================================
ROOT eigenvalues: 
------------------------------------
State configuration | energy eigenvalue (eV)
------------------------------------
     1s^{1/2} J:1/2 | -1.358351e+01
------------------------------------
     2s^{1/2} J:1/2 | -3.400454e+00
------------------------------------
     3s^{1/2} J:1/2 | -1.511596e+00
------------------------------------
     4s^{1/2} J:1/2 | -8.503202e-01
------------------------------------
     5s^{1/2} J:1/2 | -5.442173e-01
------------------------------------
     6s^{1/2} J:1/2 | -3.779333e-01
------------------------------------
     7s^{1/2} J:1/2 | -2.776697e-01
------------------------------------
     8s^{1/2} J:1/2 | -2.126026e-01
------------------------------------
     9s^{1/2} J:1/2 | -1.680162e-01
------------------------------------
    10s^{1/2} J:1/2 | -1.361489e-01
------------------------------------
    11s^{1/2} J:1/2 | -1.125290e-01
------------------------------------
    12s^{1/2} J:1/2 | -9.429455e-02
------------------------------------
    13s^{1/2} J:1/2 | -7.808953e-02
------------------------------------
    14s^{1/2} J:1/2 | -6.173130e-02
------------------------------------
    15s^{1/2} J:1/2 | -4.672159e-02
------------------------------------
    16s^{1/2} J:1/2 | -3.735581e-02
------------------------------------
    17s^{1/2} J:1/2 | -2.954559e-02
------------------------------------
    18s^{1/2} J:1/2 | -1.918991e-02
------------------------------------
    19s^{1/2} J:1/2 | -1.069129e-02
------------------------------------
    20s^{1/2} J:1/2 | -4.676752e-03
------------------------------------
    21s^{1/2} J:1/2 | -1.147953e-03
------------------------------------
=================================================
Output ROOT file: eigenvectors.root
=================================================
```

## Plot eigenvectors as a function of the radius:
```
root -l eigenvectors.root 
root [0] 
Attaching file eigenvectors.root as _file0...
root [1] graph_N8->GetXaxis()->SetRangeUser(0.0, 100.0); 
root [2] graph_N8->Draw("A*L");                          
```

