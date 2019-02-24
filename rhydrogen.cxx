#include <stdio.h>
#include <iomanip>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "rhydrogen.h"

Hydrogen::Hydrogen(Int_t N, Double_t A, Double_t s, Int_t l, Double_t j, Int_t Z, Int_t Option){
    m_A = A;
    m_s = s;
    m_N = N;
    m_Z = (Double_t)Z;
    m_Dim = 2*m_N;
    m_j = j;
    m_l = l;
    if(m_j == m_l + 0.5){
	m_Kappa = -(m_l + 1.0);
    }else if(m_j == m_l - 0.5){
	m_Kappa = m_l;
    }else{
	cerr << "ERROR: Wrong l or j !!" << endl;
	exit(-1);
    }

    if(m_N < 20){
      cerr << "ERROR: Number of grid points " << m_N << " too small! Try with 20 < N < 200" << endl;
      exit(-1);
    }

    if(m_l > 15){
      cerr << "ERROR: too high orbital quantum number " << m_l << " consult author for such feature." << endl;
      exit(-1);
    }

    m_Option = Option;
    m_Grid = new TArrayD(m_N);
    m_Ham_ar = new TArrayD(m_Dim*m_Dim);
    m_Hamilton = new TMatrixT<double>(m_Dim, m_Dim);
    m_gX = new TArrayD(m_N);

    m_ev_ar = NULL;
    m_evecs_ar = NULL;

    state_orb_map[0.0] = "s";
    state_orb_map[1.0] = "p";
    state_orb_map[2.0] = "d";
    state_orb_map[3.0] = "f";
    state_orb_map[4.0] = "g";
    state_orb_map[5.0] = "h";
    state_orb_map[6.0] = "l";
    state_orb_map[7.0] = "m";
    state_orb_map[8.0] = "n";
    state_orb_map[9.0] = "o";
    state_orb_map[10.0] = "p";
    state_orb_map[11.0] = "q";
    state_orb_map[12.0] = "r";
    state_orb_map[13.0] = "s";
    state_orb_map[14.0] = "t";
    
    state_spin_map[0.5] = "1/2";
    state_spin_map[1.5] = "3/2";
    state_spin_map[2.5] = "5/2";
    state_spin_map[3.5] = "7/2";
    state_spin_map[4.5] = "9/2";
    state_spin_map[5.5] = "11/2";

    cout << "State: n " << state_orb_map[(double)m_l] << "^{" << state_spin_map[TMath::Abs((double)(m_j - m_l))] << "}" << endl;

}

Hydrogen::~Hydrogen(){

}

Double_t Hydrogen::R_theta(Int_t i){
//r(theta) = [s*theta - A*arctan(s*theta/A)]/(Pi - theta)^{2}
//
//Remember: theta == m_Grid->At(i) since we are on a grid

    if(i >= m_N)cerr << "In R_theta(i) i is greater than N - this is not OK" << endl;

    return (((m_s * m_Grid->At(i)) - (m_A*atan(m_s*m_Grid->At(i)/m_A)))/((PI - m_Grid->At(i))*(PI - m_Grid->At(i))));
}

Double_t Hydrogen::Jacobian(Int_t i){
//J(theta) = dr/dtheta = 
//    [s - (s/(1+((s*theta)/A)^{2}))]/[(Pi - theta)^{2}] +
//     +  2 * [r(theta)]/[(Pi-theta)]
    
    if(i >= m_N)cerr << "In Jacobian(i) i is greater than N - this is not OK" << endl;

    return ((m_s - (m_s/(1.0+TMath::Power(m_s*m_Grid->At(i)/m_A, 2.0))))/(TMath::Power(PI - m_Grid->At(i),2.0)) + (2.0*R_theta(i)/(PI - m_Grid->At(i))));
}

Double_t Hydrogen::CosineSineMatrix(Int_t i, Int_t j){

    if(i >= m_N || j >= m_N)cerr << "In CosineSineMatrix(i,j) i or j is greater than N - this is not OK" << endl;
    Double_t ret = 0.0;
    for(Int_t m = 1; m <= m_N; m++){
	ret = ret + (m*cos(m*m_Grid->At(i))*sin(m*m_Grid->At(j)));
    }
   ret *= (2.0/(m_N+1.0));
    return ret;
}


Double_t Hydrogen::Delta(Int_t i, Int_t j){
    Double_t ret = 0.0;
    if(i == j)ret = 1.0;
    return ret;
}

//Setup functions
void Hydrogen::SetupThetaGrid(){
// \delta theta = Pi/N+1

    for(Int_t i = 1; i <= m_N; i++){
	m_Grid->SetAt((i*PI)/(m_N+1.0), i-1);
	m_gX->SetAt(m_Grid->At(i-1), i-1);
//	cout << m_Grid->At(i) << endl;
    }
}



void Hydrogen::CalculateHamiltonian(){
//Calculates the Hamiltonian

    Double_t t_Kappa = 0.0;
    for(int i = 0; i < m_Dim; i++){
	for(int j = 0; j < m_Dim; j++){

	    if(i < m_N && j < m_N){//upper left block - diagonal
	      m_Ham_ar->SetAt( -factor()*((m_Z*ALPHA/R_theta(i)) - 1.0)*Delta(i,j) , j+i*m_Dim);

	    }else if(i < m_N && j >= m_N){//upper right block - off-diagonal
		m_Ham_ar->SetAt( -factor()*((CosineSineMatrix(i,j-m_N)/Jacobian(i)) - (m_Kappa*Delta(i,j-m_N)/R_theta(i))), j+i*m_Dim);

	    }else if(i >= m_N && j < m_N){//lower left block - off-diagonal
	      m_Ham_ar->SetAt( factor()*((CosineSineMatrix(i-m_N,j)/Jacobian(i-m_N)) + (m_Kappa*Delta(i-m_N,j)/R_theta(i-m_N))), j+i*m_Dim);

	    }else if(i >= m_N && j >= m_N){//lower right block - diagonal
	      m_Ham_ar->SetAt( -factor()*((m_Z*ALPHA/R_theta(i-m_N)) + 1.0)*Delta(i-m_N, j-m_N) , j+i*m_Dim);
	    }
	}
    }

   m_Hamilton->SetMatrixArray(m_Ham_ar->GetArray());

//	(*m_Hamilton) *= -1.0;//maybe ...
//    m_Hamilton->Print();

}


void Hydrogen::DiagonalizeHamiltonian(){
    
  TMatrixDEigen m_in(*m_Hamilton);
  TMatrixD m_ev = m_in.GetEigenValues();
//    m_ev.Print();
  TMatrixD m_evecs = m_in.GetEigenVectors();
//    m_evecs.Print();

    this->PrintEigenValues(&m_ev);

    Double_t * m_ev_ar = m_ev.GetMatrixArray();
    Double_t * m_evecs_ar = m_evecs.GetMatrixArray();

    this->SaveEigenVectors(m_ev_ar, m_evecs_ar);

    this->TestUnityOrth(m_ev_ar, m_evecs_ar);


    /*
//GSL

    cout << endl;
    cout << "GSL eigenvalues: " << endl;
    cout << "_________________" << endl;


    gsl_matrix_view m = gsl_matrix_view_array (m_Ham_ar->GetArray(), m_Dim, m_Dim);

//print the results    
//     gsl_matrix * mm = (gsl_matrix*) &m;
//     for (int i = 0; i < m_Dim; i++){
// 	for (int j = 0; j < m_Dim; j++){
// 	    printf ("%g   ", gsl_matrix_get (mm, i, j));
// 	}
// 	printf("\n");
//     }


    gsl_vector_complex *eval = gsl_vector_complex_alloc (m_Dim);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (m_Dim, m_Dim);
     
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (m_Dim);
    
    gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
    
    gsl_eigen_nonsymmv_free (w);
    
    gsl_eigen_nonsymmv_sort (eval, evec, 
			  GSL_EIGEN_SORT_ABS_ASC);
       
    {
	int i;
	

	for (i = 0; i < m_Dim; i++)
	{
	    gsl_complex eval_i 
                = gsl_vector_complex_get (eval, i);
	    gsl_vector_complex_view evec_i 
                = gsl_matrix_complex_column (evec, i);

//		    printf ("eigenvalue = %g + %gi\n",
//			    GSL_REAL(eval_i), GSL_IMAG(eval_i));
	    if(TMath::Abs(GSL_REAL(eval_i)) < 1.0){
		if(GSL_REAL(eval_i) > 0.0){
		    printf ("eigenvalue = %1.12E + %Ei  eV \n",
			    (GSL_REAL(eval_i)*EMASS_EV - EMASS_EV), GSL_IMAG(eval_i));
		}
	    }
//	    if(i == 0){
//		m_base = eval_i - EMASS*SOL*SOL/EV;
//		printf ("eigenvalue = %g ---> 0.000 eV (n=%d) \n", m_base, i+1 );
//	    }
//	    if(i > 0){
//		printf ("eigenvalue = %g ---> %g eV (n=%d)\n", eval_i - EMASS*SOL*SOL/EV, (eval_i - (EMASS*SOL*SOL/EV)- m_base), i+1);
//	    }
//	    printf ("eigenvector = \n");
//	    gsl_vector_fprintf (stdout, 
//				&evec_i.vector, "%g");
	}
    }
    
    gsl_vector_complex_free (eval);
    gsl_matrix_complex_free (evec);
    */
    
}

void Hydrogen::PrintEigenValues(TMatrixD * ev){

  //eigenvalues
  printf("Electron mass: %E eV\n", EMASS_EV);

    cout << endl;
    cout << "=================================================" << endl;
    cout << "ROOT eigenvalues of bound states: " << endl;

    Double_t * m_ev_ar = ev->GetMatrixArray();
    map<int, double> qn_evals;
    int count_boundstates = 0;
    cout << setw(10) << "------------------------------------" << endl;
    cout << setw(10) << "State configuration" << setw(3) << " | " << setw(20) << "energy eigenvalue (eV)" << endl;
    cout << setw(10) << "------------------------------------" << endl;
    int PQN = 0;
    for(int i = m_Dim-1; i >= 0; i--){
 	if(TMath::Abs(m_ev_ar[i+i*m_Dim]) < 1.0){
 	    if(m_ev_ar[i+i*m_Dim] > 0.0){
	      qn_evals[count_boundstates] = m_ev_ar[i+i*m_Dim]*EMASS_EV - EMASS_EV;
	      count_boundstates++;
	    }
	}
    }

    int PrincQN = 0;
    Double_t evalue = 0.0;
    const int size = qn_evals.size();
    for(int count = 0; count < size; count++){
      //    for(map<int, double>::iterator it = qn_evals.begin(); it != qn_evals.end(); it++){
      PrincQN = (int)m_l + 1 + count;
      evalue = qn_evals[count];
      if(m_j < m_l){
	evalue = qn_evals[count+1];
	if(count == size -1){
	  break;
	}
      }
      cout << setw(6) << PrincQN << state_orb_map[(double)m_l] << "^{" << state_spin_map[TMath::Abs((double)(m_j - m_l))] << "} J:" << state_spin_map[m_j] << " | " << setw(10) << scientific << evalue << endl;
      cout << setw(10) << "------------------------------------" << endl;
    }

}


void Hydrogen::SaveEigenVectors(Double_t * m_ev_ar, Double_t * m_evecs_ar ){

    TFile * file = new TFile("eigenvectors.root", "recreate");

    const int n_arrs = 10;
    TGraph ** gr = new TGraph*[n_arrs];
    TGraph ** gr2 = new TGraph*[n_arrs];
    Double_t ** arr = new Double_t*[n_arrs];
    Double_t ** arr2 = new Double_t*[n_arrs];
    Double_t * grid = new Double_t[(m_Dim/2)-1];
    for(int i = 0; i < n_arrs; i++){
	arr[i] = new Double_t[(m_Dim/2)-1];
	arr2[i] = new Double_t[(m_Dim/2)-1];
    }
    
    int count_boundstates = -1;
    for(int i = m_Dim-1; i >= 0; i--){
	if(count_boundstates == n_arrs-1)
	    break;
	if(TMath::Abs(m_ev_ar[i+i*m_Dim]) < 1.0){
	    if(m_ev_ar[i+i*m_Dim] > 0.0){
		count_boundstates++;
//		cout << "boundstates: " << count_boundstates << endl;
		//read columns of m_evecs matrix
		for(int r = (m_Dim/2)-1-1; r >= 0; r--){
		    arr2[count_boundstates][r] = TMath::Power(m_evecs_ar[i+r*m_Dim], 2.0);
		    arr[count_boundstates][r] = m_evecs_ar[i+r*m_Dim];
		    grid[r] = R_theta(r);
		    //		    cout << grid[r] << ", " << arr2[count_boundstates][r] << endl;
		}
		//	cout << endl;
	    }
	}
    }


    string prefix = "graph";
    string name = "";
    int PrincQN = 0;
    int j = 0;
    for(int i = 0; i < n_arrs; i++){
	PrincQN = (int)m_l + 1 + i;
	if(m_j < m_l && i < n_arrs-1)
	  j = i+1;
	else if(m_j < m_l && i == n_arrs-1)
	  break;
	else if(m_j > m_l)
	  j = i;
	gr2[i] = new TGraph((m_Dim/2)-1, grid, arr2[j]);
	gr[i] = new TGraph((m_Dim/2)-1, grid, arr[j]);
	name = prefix + "_N" + string(itoa(PrincQN));
	gr[i]->SetName(name.c_str());
	name = prefix + "2_N" + string(itoa(PrincQN));
	gr2[i]->SetName(name.c_str());
	gr[i]->Write();
	gr2[i]->Write();
    }

    file->Close();
    delete file;
    delete[] arr;
    delete[] arr2;
    delete[] gr;
    delete[] gr2;
    delete grid;
    cout << "=================================================" << endl;
    cout << "Output ROOT file: eigenvectors.root" << endl;



}


 void Hydrogen::TestUnityOrth(Double_t * m_ev_ar, Double_t * m_evecs_ar){

    const int n_arrs = 10;
    Double_t ** arr = new Double_t*[n_arrs];
    Double_t ** arr2 = new Double_t*[n_arrs];
    for(int i = 0; i < n_arrs; i++){
	arr[i] = new Double_t[(m_Dim/2)-1];
	arr2[i] = new Double_t[(m_Dim/2)-1];
    }


    int count_boundstates = -1;
    for(int i = m_Dim-1; i >= 0; i--){
	if(count_boundstates == n_arrs-1)
	    break;
	if(TMath::Abs(m_ev_ar[i+i*m_Dim]) < 1.0){
	    if(m_ev_ar[i+i*m_Dim] > 0.0){
		count_boundstates++;
//		cout << "boundstates: " << count_boundstates << endl;
		//read columns of m_evecs matrix
		for(int r = (m_Dim/2)-1-1; r >= 0; r--){
		    arr2[count_boundstates][r] = TMath::Power(m_evecs_ar[i+r*m_Dim], 2.0);
		    arr[count_boundstates][r] = m_evecs_ar[i+r*m_Dim];
		    //		    cout << grid[r] << ", " << arr2[count_boundstates][r] << endl;
		}
		//	cout << endl;
	    }
	}
    }




    //orthogonality?
    cout << "=================================================" << endl;
    cout << "Test orthogonality and unity for first 2 vectors:" << endl;

    Double_t proj = 0.0;
    Double_t a1 = 0.0;
    Double_t a2 = 0.0;
    for(int i = 0; i < (m_Dim/2)-1; i++){
	proj += arr[0][i]*arr[1][i];
	a1 += arr[0][i]*arr[0][i];
	a2 += arr[1][i]*arr[1][i];
    }
    a1 = TMath::Sqrt(a1);
    cout << "_______________________" << endl;
    cout << "Magnitude of vectors: " << endl;
    cout << "vector1: " << a1 << endl;
    a2 = TMath::Sqrt(a2);
    cout << "vector2: " << a2 << endl;
    proj = TMath::ACos(proj/(a1*a2));    
    //180 = Pi radian --> 1 = Pi/180 radian, 1 radian = 180/Pi 
    proj = proj*180.0/TMath::Pi();
    cout << "_______________________" << endl;
    cout << "Projection vector1 on vector2: " << proj << " degrees" <<  endl;
    cout << "=================================================" << endl;

    delete[] arr;
    delete[] arr2;

 }


Double_t Hydrogen::factor(){

//    return EMASS*SOL*SOL/EV;//is this okay???

    return 1.0;//is this okay???
}

//Utility functions

char * Hydrogen::itoa(int x){
    char * res = new char[128];
    sprintf(res, "%d", x);
    return res;
}



//Note: in ROOT to fill up a TMatrix column wise one must
//do e.g.:
//
// TArrayD * cm = new TArrayD(25); //5x5 matrix
// TMatrixT<double> * m55 = new TMatrixT<double>(5,5);
// Double_t temp = 0.0;
// for(int i = 0; i < 5; i++){
//   for(int j = 0; j < 5; j++){
//      cm->SetAt(temp, i+j*5);
//      temp += 1.0;
//   }
// }
// m55->SetMatrixArray(cm->GetArray());
// m55->Print();
// 5x5 matrix is as follows
//
//    |      0    |      1    |      2    |      3    |      4    |
//----------------------------------------------------------------------
//   0 |          0           5          10          15          20 
//   1 |          1           6          11          16          21 
//   2 |          2           7          12          17          22 
//   3 |          3           8          13          18          23 
//   4 |          4           9          14          19          24 
//
//don't ask why I've put this into here....
//
