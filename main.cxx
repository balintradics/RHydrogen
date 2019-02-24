#include "rhydrogen.h"
#include "time.h"
#include "cmdline.h"
#include <iomanip>

int main(int argc, char * argv[]){

  if(argc == 1){
    cout << "------------------------------------------------------------------" << endl;
    cout << "|                     R H Y D R O G E N                          |" << endl;
    cout << "| RHydrogen: this c++ program solves the Dirac equation for      |" << endl ;
    cout << "| simple Hydrogen-like atoms, and calculates the eigenvectors    |" << endl ;
    cout << "| and eigenvalues of orbitals given the quantum numbers of the   |" << endl;
    cout << "| particular state.                                              |" << endl;
    cout << "|                                                                |" << endl;
    cout << "| It is based on a general numerical differential equation       |" << endl;
    cout << "| solving algorithm inspired by the paper:                       |" << endl;
    cout << "| E.Ackad and M. Horbatsch: Numerical solution of the            |" << endl;
    cout << "| Dirac equation by a mapped Fourier grid method                 |" << endl;
    cout << "| J.Phys.A:Math.Gen. 38 (2005) 3157.                             |" << endl;
    cout << "|                                                                |" << endl;
    cout << "| Program author: B. Radics (Balint.Radics@cern.ch)              |" << endl;
    cout << "|                                                                |" << endl;
    cout << "| Try to launch with some options                                |" << endl ;
    cout << "| type: rhydrogen --help for the complete list of options        |" << endl ;
    cout << "------------------------------------------------------------------" << endl;
    exit(-1);
  }

  TStopwatch timer;
  timer.Start();


  gengetopt_args_info args_info;
  
  Int_t l;
  Double_t j;
  Int_t Z;
  Double_t A;
  Double_t s;
  Int_t N;
  Int_t Op = 0; //0 -- F(x) = 1, 1 -- F(x) = something..not working yet
 

  /* let's call our cmdline parser */
  cmdline_parser (argc, argv, &args_info);
  //  if (cmdline_parser (argc, argv, &args_info) != 0)
  //  exit(-1) ;

  if (!args_info.number_gridpts_given || !args_info.scale_param_given || !args_info.l_orbital_qn_given || !args_info.j_total_qn_given || !args_info.Z_nucl_charge_given){
    cerr << "Error: no parameters given!" << endl;
    exit(-1);
  }else{
    s = args_info.scale_param_arg;
    N = args_info.number_gridpts_arg;
    l = args_info.l_orbital_qn_arg;
    j = args_info.j_total_qn_arg;
    Z = args_info.Z_nucl_charge_arg;
  }

  cout << "Scale parameter: " << s << endl;
  cout << "Grid points: " << N << endl;
  cout << "Orbital angular mom. quantum number: " << l << endl;
  cout << "Total angular mom. quantum number: " << j << endl;
  cout << "Z nucl. charge: " << Z << endl;

  cmdline_parser_free (&args_info); /* release allocated memory */


    Hydrogen * df = new Hydrogen(N, A, s, l, j, Z, Op);
    df->SetupThetaGrid();
    df->CalculateHamiltonian();
    df->DiagonalizeHamiltonian();

    delete df;
//     Int_t Ni = 300;
//     TArrayD * s_arr_20 = new TArrayD(Ni);
//     TArrayD * e_arr_20 = new TArrayD(Ni);
//     TArrayD * s_arr_40 = new TArrayD(Ni);
//     TArrayD * e_arr_40 = new TArrayD(Ni);
//     TArrayD * s_arr_60 = new TArrayD(Ni);
//     TArrayD * e_arr_60 = new TArrayD(Ni);

//     Hydrogen * df;
//     Double_t m_s = 0;


//     for(int i = 0; i < Ni; i++){
// 	m_s = 100.0 + i*300.0;
// 	df = new Hydrogen(100, A, m_s, l, j, Z, Op);
// 	df->SetupThetaGrid();
// 	df->CalculateHamiltonian();
// 	df->DiagonalizeHamiltonian();
// 	s_arr_20->SetAt(m_s, i);
// 	e_arr_20->SetAt(df->GetNRG(N), i);
// 	delete df;
// 	df = new Hydrogen(200, A, m_s, l, j, Z, Op);
// 	df->SetupThetaGrid();
// 	df->CalculateHamiltonian();
// 	df->DiagonalizeHamiltonian();
// 	s_arr_40->SetAt(m_s, i);
// 	e_arr_40->SetAt(df->GetNRG(N), i);
// 	delete df;
// 	df = new Hydrogen(300, A, m_s, l, j, Z, Op);
// 	df->SetupThetaGrid();
// 	df->CalculateHamiltonian();
// 	df->DiagonalizeHamiltonian();
// 	s_arr_60->SetAt(m_s, i);
// 	e_arr_60->SetAt(df->GetNRG(N), i);
// 	if(i == Ni - 1)cout << df->GetNRG(N) << endl;
// 	delete df;



//     }

//     TFile * file = new TFile("out.root", "recreate");
    
//     TGraph * g_se_20 = new TGraph(Ni, s_arr_20->GetArray(), e_arr_20->GetArray());
//     g_se_20->SetNameTitle("N_20", "N_20");
//     TGraph * g_se_40 = new TGraph(Ni, s_arr_40->GetArray(), e_arr_40->GetArray());
//     g_se_40->SetNameTitle("N_40", "N_40");
//     TGraph * g_se_60 = new TGraph(Ni, s_arr_60->GetArray(), e_arr_60->GetArray());
//     g_se_60->SetNameTitle("N_60", "N_60");

//     g_se_20->Write();
//     g_se_40->Write();
//     g_se_60->Write();
//     delete file;

    timer.Stop();
    cout << "Real time: " << setw(5) << timer.RealTime() << " s" << endl; 
    cout << "CPU time:  " << setw(5) << timer.CpuTime() << " s" << endl;
    cout << "=================================================" << endl;
    cout << "To check reality of results, please consult the NIST Atomic Database: " << endl;
    cout << endl;
    cout << "    http://physics.nist.gov/PhysRefData/ASD/levels_form.html" << endl;
    cout << endl;
    cout << "=================================================" << endl;
    if(Z > 20){
      cout << "\nWARNING: You specified Z = " << Z << " but accuracy of calculation decreases with increasing Z," << endl;
      cout << "because nuclear-electron spin/magnetic/etc. interactions are not yet included...\n" << endl;
    }
    cout << "Bye-bye" << endl;
    cout << endl;

    return 0;

 }


