#include "TMatrixD.h"
#include "TRandom3.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TDecompSVD.h"
#include "TDecompChol.h"
#include "TLine.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TPad.h"

const int nthrows = 1;

class MatrixCollect {

  public:
    MatrixCollect(TMatrixDSym *input);
    ~MatrixCollect() {};
    void ConstructNew();
    void Write(std::string);
    // How many eigen values for a cut
    int HowMany(double cut);
    void SetLimits(TH2D*);
    void TestMatrix();

  private:
    // The eigenmatrix
    TMatrixDSymEigen *eigen;
    // The original matrix
    TMatrixDSym *origmat;
    // The eigenvalues
    //TVectorD *eigenval;
    // The eigenvectors
    //TMatrixD *eigenvec;
    // Number of parameters
    int npars;

    int nvals;

    // Number of cuts
    int ncuts;
    // Vector of cuts
    std::vector<int> cuts;
    // The transfer matrix
    std::vector<TMatrixD> transfer;
    // The transposed transfer matrix
    std::vector<TMatrixD> transferT;
    // The updated matrix
    std::vector<TMatrixD*> newmat;

};

// Clone the matrices, set up eigen matrix
MatrixCollect::MatrixCollect(TMatrixDSym *input) {

  std::cout << "Input matrix with " << input->GetNrows() << "x" << input->GetNcols() << " (rows x columns)" << std::endl;
  //input->Print();

  origmat = new TMatrixDSym(*input);
  //origmat.SetMatrixArray(input->GetMatrixArray());

  // Then construct the TMatrixDEigen
  eigen = new TMatrixDSymEigen(*input);

  // Get the eigen values
  //std::cout << eigen->GetEigenValues() << std::endl;

  //eigenval = &(eigen->GetEigenValues());
  //eigenvec = &(eigen->GetEigenVectors());

  // The number of parameters
  npars = origmat->GetNrows();

  std::cout << "Done constructing" << std::endl;
}

int MatrixCollect::HowMany(double cut) {
  // Now loop through looking for small eigen values
  int nvals_eigen = 0;
  for (int i = 0; i < npars; ++i) {
    if (eigen->GetEigenValues()(i) > cut) nvals_eigen++;
    if (nvals_eigen > npars) break;
  }
  return nvals_eigen;
}

void MatrixCollect::SetLimits(TH2D *plot) {
  double maximum = plot->GetMaximum() > fabs(plot->GetMinimum()) ? plot->GetMaximum() : fabs(plot->GetMinimum());
  plot->SetMaximum(maximum);
  plot->SetMinimum(-1*maximum);
}

// Construct the reduced matrix
void MatrixCollect::ConstructNew() {

  // Convenience
  TMatrixD EigenVectors = eigen->GetEigenVectors();
  TVectorD EigenValues = eigen->GetEigenValues();

  ncuts = 5;

  // Loop over the cuts
  for (int i = 0; i < ncuts; ++i) {

    nvals = (i+1)*npars/ncuts;
    std::cout << "Cutting on " << nvals << " eigen values" << std::endl;

    // The reduced eigenvector matrix
    TMatrixD transfer_temp(EigenVectors.GetSub(0, EigenVectors.GetNcols()-1, 0, nvals-1));
    transfer.push_back(transfer_temp);

    // The transposed transfer matrix
    TMatrixD transferT_temp = transfer_temp;
    transferT_temp.T();
    transferT.push_back(transferT_temp);

    // The eigen value matrix
    TMatrixD reduced(nvals, nvals);
    reduced.Zero();
    for (int i = 0; i < nvals; ++i) {
      reduced(i,i) = EigenValues(i);
    }

    // The recreated matrix
    TMatrixD temp = transfer_temp*reduced;
    TMatrixD temp2 = temp*transferT_temp;

    TMatrixD *newmat_temp = new TMatrixD(temp2);

    // push back the final matrix
    newmat.push_back(newmat_temp);
  } // Finished looping over the cuts

  std::cout << "Finished constructing " << ncuts << " new covariances" << std::endl;
}

// Try to build the uncertainty band
void MatrixCollect::TestMatrix() {

  // Random number generator
  TRandom3 *RandNo = new TRandom3(1337);

  // The vector of random numbers
  TVectorD *randoms = new TVectorD(npars);

  // First make a Cholesky decomp of the original matrix
  TDecompChol chol(*origmat);
  if (!chol.Decompose()) {
    std::cout << "Failed to decompose new matrix" << std::endl;
  } else {
    std::cout << "Decomposed new matrix" << std::endl;
  }
  TMatrixD cholmat(chol.GetU());

  // Propagate the throw through the transfer matrix to the eigenbasis
  TMatrixD main_transfer = transfer.back();

  // Loop over the throws and make 
  for (int i = 0; i < nthrows; ++i) {
    //if (i%(nthrows/5) == 0) std::cout << "On throw " << i << "/" << nthrows << std::endl;
    // Make Gaussian variations around the nominal
    for (int j = 0; j < npars; ++j) {
      (*randoms)(j) = RandNo->Gaus(0, 1);
    }
    // The throw of parameters
    TVectorD newthrow = cholmat*(*randoms);

    // The throw of the parameters in the eigen basis
    TVectorD newthrow_basis = main_transfer*newthrow;

    newthrow_basis.Print();

    // Now have a throw in the eigenbasis
  }

}

// Want to write
// Original matrix
// Eigen-values
// For each cut
//  - recreated original matrix
//  - throws
void MatrixCollect::Write(std::string filename) {

  // The output file
  TFile *output = new TFile(Form("%s_test.root", filename.c_str()), "recreate");

  // Make a TH2D of the original matrix
  TH2D *OrigTH2 = new TH2D("Original_Matrix", "Original_Matrix", npars, 0, npars, npars, 0, npars); 
  for (int i = 0; i < npars; ++i) {
    for (int j = 0; j < npars; ++j) {
      OrigTH2->SetBinContent(i+1, j+1, (*origmat)(i, j));
    }
  }
  OrigTH2->Write("Original_Matrix");

  // Write a TH1D of the eigenvalues
  TH1D *EigenVals = new TH1D("Eigen_Values", "Eigen_Values", npars, 0, npars);
  for (int i = 0; i < npars; ++i) {
    EigenVals->SetBinContent(i+1, (eigen->GetEigenValues())(i));
  }
  EigenVals->Write("Eigen_Values");

  // Loop over each of the new matrices
  TH2D *temp_th2 = new TH2D("temp_th2", "temp_th2", npars, 0, npars, npars, 0, npars);
  for (int i = 0; i < ncuts; ++i) {
    TMatrixD *temp_mat = newmat.at(i);
    temp_th2->Reset();
    temp_th2->SetMinimum(-1111);
    temp_th2->SetMaximum(-1111);
    TH2D *temp_th2_copy = (TH2D*)temp_th2->Clone();
    temp_th2_copy->SetMinimum(-1111);
    temp_th2_copy->SetMaximum(-1111);
    for (int j = 0; j < npars; ++j) {
      for (int k = 0; k < npars; ++k) {
        temp_th2->SetBinContent(j+1, k+1, (*temp_mat)(j,k));
        temp_th2_copy->SetBinContent(j+1, k+1, (temp_th2->GetBinContent(j+1, k+1)-OrigTH2->GetBinContent(j+1, k+1))/OrigTH2->GetBinContent(j+1, k+1));
      }
    }

    temp_th2->Write(Form("Recreated_%i", i));
    // Copy
    temp_th2->Divide(OrigTH2);
    temp_th2->SetMinimum(0.5);
    temp_th2->SetMaximum(1.5);
    temp_th2_copy->SetMinimum(-0.1);
    temp_th2_copy->SetMaximum(0.1);
    temp_th2->Write(Form("Recreated_%i_RelNom", i));
    // Make percentage difference
    temp_th2_copy->Write(Form("Recreated_%i_RelNomPercent", i));
  }

  // Close the file
  std::cout << "Wrote output to " << output->GetName() << std::endl;
  output->Close();
}

void makedecomp(std::string filename) {

  // The input file
  TFile *file = new TFile(filename.c_str());

  // The collection of matrices
  TMatrixDSym *mat = (TMatrixDSym*)file->Get("Covariance_Matrix")->Clone();

  // Create the large object
  MatrixCollect matobj(mat);

  // Construct all the new matrices
  matobj.ConstructNew();

  matobj.TestMatrix();

  // Write to file
  matobj.Write(filename);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cout << "Need one argument, the filename" << std::endl;
    return -1;
  }

  makedecomp(std::string(argv[1]));
}
