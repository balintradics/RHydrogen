/////////////////////////////////////////////////////////////////////
//  General numeric differential equation solver class
//  base on the paper:
//  E.Ackad and M. Horbatsch: Numerical solution of the Dirac equation
//              by a mapped Fourier grid method
//  J.Phys.A:Math.Gen. 38 (2005) 3157.
/////////////////////////////////////////////////////////////////////

#ifndef HYDROGEN_H
#define HYDROGEN_H

#include<TArrayD.h>
#include<TMatrixD.h>
#include<TMatrixDEigen.h>
#include<TVectorD.h>
#include<TFile.h>
#include<TF1.h>
#include<TGraph.h>
#include<TStopwatch.h>
#include<math.h>
#include<TMath.h>
#include<iostream>
#include<fstream>
#include<string>
#include<stdlib.h>
#include<map>

using namespace std;

//Physical Constants
#define HBAR 1.054571628e-34 //J*s
#define ALPHA 7.2973525376e-3
#define SOL 299792458 //m*s^{-1}
#define EV 1.602176487e-19 //J
#define EMASS 9.10938215e-31 //kg
#define EMASS_EV EMASS*SOL*SOL/EV //eV
#define PERM 8.854187817e-12 //F*m^{-1}
#define PI 3.14159265358979312
#define HARTREE 4.35974394e-18 //J
#define HARTREE_EV 27.21138386 //eV
#define KG_HARTREE 2.06148616e+34 //Eh

class Hydrogen{

 public:
    Hydrogen(Int_t N, Double_t A, Double_t s, Int_t l, Double_t j, Int_t Z, Int_t Option);
    virtual ~Hydrogen();

//Setup functions
    Double_t R_theta(Int_t i);
    Double_t Jacobian(Int_t i);
    Double_t Delta(Int_t i, Int_t j);
    Double_t CosineSineMatrix(Int_t i, Int_t j);
    void SetupThetaGrid();
    void CalculateHamiltonian();//Calculates the Hamiltonian
    void DiagonalizeHamiltonian();//Diagonalizes the Hamiltonian
    void PrintEigenValues(TMatrixD * ev);
    void SaveEigenVectors(Double_t * m_ev_ar, Double_t * m_evecs_ar);
    void TestUnityOrth(Double_t * m_ev_ar, Double_t * m_evecs_ar);
    Double_t GetNRG(Int_t n){return m_base;}

//utility functions
    char * itoa(int x);
    Double_t factor();//unit conversion

 public:
    TArrayD * m_Grid;//the discrete playground in theta
    TArrayD * m_Ham_ar;
    TMatrixT<double> * m_Hamilton;//Hamiltonian matrix
    TArrayD * m_gX;//array for graph, unit: hbar/mc
    TArrayD * m_gY;//array for graph
    Double_t m_s;// some scale variable if needed
    Double_t m_A;// some scale variable if needed
    Double_t * m_ev_ar;
    Double_t * m_evecs_ar;
    Int_t m_N; // grid points
    Double_t m_Kappa;// -(l+1) for j = l+1/2; l for j = l - 1/2
    Double_t m_j;
    Double_t m_base;
    Int_t m_l;
    Double_t m_Z;
    Int_t m_Dim;//dimension of matrices
    Int_t m_Option;//Option on F(x)
    map<double,string> state_orb_map;
    map<double,string> state_tot_map;
    map<double,string> state_spin_map;

};

#endif
