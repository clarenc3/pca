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

const int nthrows = 20000;

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
    TVectorD *eigenval;
    // The eigenvectors
    TMatrixD *eigenvec;

    // Number of parameters
    int npars;

    // Number of cuts
    int ncuts;

    // The transfer matrices for all the cuts
    std::vector<TMatrixD> transfer;

    // The transposed transfer matrices for all the cuts
    std::vector<TMatrixD> transferT;

    // The updated matrix for all the cuts
    std::vector<TMatrixD*> newmat;

    // The histograms for each of the cuts
    std::vector<TH2D*> Throws;
};

// Clone the matrices, set up eigen matrix
MatrixCollect::MatrixCollect(TMatrixDSym *input) {

  std::cout << "Input matrix with " << input->GetNrows() << "x" << input->GetNcols() << " (rows x columns)" << std::endl;

  // Clone the matrix
  origmat = new TMatrixDSym(*input);

  // Then construct the TMatrixDEigen
  eigen = new TMatrixDSymEigen(*input);

  // Copy over the eigen values and vectors
  eigenval = new TVectorD(eigen->GetEigenValues());
  eigenvec = new TMatrixD(eigen->GetEigenVectors());

  // The number of parameters in total
  npars = origmat->GetNrows();

  std::cout << "Done constructing eigen matrix, values and vectors" << std::endl;
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

  // Put 5 cuts in
  ncuts = 10;

  // Loop over the cuts
  for (int i = 0; i < ncuts; ++i) {
    // Number of eigen values for this cut
    int nvals = (i+1)*npars/ncuts;
    std::cout << "Cutting on " << nvals << " eigen values" << std::endl;

    // The reduced eigenvector matrix
    TMatrixD transfer_temp(eigenvec->GetSub(0, eigenvec->GetNcols()-1, 0, nvals-1));

    //std::cout << "Transfer:" << std::endl;
    //transfer_temp.Print();
    //std::cout << "eigenvec:" << std::endl;
    //eigenvec->Print();
    //throw;

    transfer.push_back(transfer_temp);

    // The transposed transfer matrix
    TMatrixD transferT_temp = transfer_temp;
    transferT_temp.T();
    transferT.push_back(transferT_temp);

    // The eigen value matrix
    TMatrixD reduced(nvals, nvals);
    reduced.Zero();
    // Write in the diagonal eigen-values 
    for (int i = 0; i < nvals; ++i) {
      reduced(i,i) = (*eigenval)(i);
    }

    // The recreated matrix
    TMatrixD *newmat_temp = new TMatrixD(transfer_temp*reduced*transferT_temp);

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
  // Get the TMatrixD of the Cholesky decomposition
  TMatrixD cholmat(chol.GetU());
  // And transpose it
  cholmat.T();

  for (int i = 0; i < ncuts; ++i) {
    TH2D *eigenthrows = new TH2D(Form("EigenThrow_%i", transfer[i].GetNcols()), Form("EigenThrow_%i;Parameter number; Parameter value", transfer[i].GetNcols()), npars, 0, npars, 1000, -2.5, 2.5);
    Throws.push_back(eigenthrows);
  }

  // Loop over the throws and make 
  for (int i = 0; i < nthrows; ++i) {

    // Be kind and verbose
    if (i%(nthrows/20) == 0) std::cout << "On throw " << i << "/" << nthrows << " (" << double(i)/nthrows*100 << "%)" << std::endl;

    // Make Gaussian variations around the nominal
    for (int j = 0; j < npars; ++j) {
      (*randoms)(j) = RandNo->Gaus(0, 1);
    }

    // The throw of parameters in the full matrix
    TVectorD newthrow = cholmat*(*randoms);

    // Now go through the eigen value cuts
    for (int j = 0; j < ncuts; ++j) {
      // The dimensions of the transfer matrix is the number of eigenvalues
      int neigen = transfer[j].GetNcols();
      // Make a smaller array of random numbers
      TVectorD randoms_eigen(neigen);
      for (int k = 0; k < neigen; ++k) randoms_eigen(k) = (*randoms)(k);
      // Now multiply this random number vector from the eigen basis to parameter basis via transfer matrix
      // And now the random taking the correlation into account
      TVectorD rands_chol_eigen = (cholmat)*(transfer[j]*randoms_eigen);
      // Loop over the parameters and fill
      for (int k = 0; k < npars; ++k) Throws[j]->Fill(k, rands_chol_eigen(k));
    } // End the cut loop
  } // End the throw loop

  // Normalise the throws
  //for (int i = 0; i < ncuts; ++i) Throws[i]->Scale(1/nthrows);

  // Delete allocated memory
  delete RandNo;
  delete randoms;
} // End the function

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
      OrigTH2->SetBinError(i+1, j+1, 0);
    }
  }
  OrigTH2->Write("Original_Matrix");

  // Let's make a diagonal only also
  TH1D *err = new TH1D("Uncertainty", "Uncertainty", npars, 0, npars);
  for (int i = 0; i < npars; ++i) err->SetBinContent(i+1, sqrt((*origmat)(i,i)));
  err->Write();

  // Write a TH1D of the eigenvalues
  TH1D *EigenVals = new TH1D("Eigen_Values", "Eigen_Values", npars, 0, npars);
  double eigentot = 0;
  for (int i = 0; i < npars; ++i) {
    EigenVals->SetBinContent(i+1, (eigen->GetEigenValues())(i));
    eigentot += (eigen->GetEigenValues())(i);
  }
  EigenVals->Write("Eigen_Values");

  // Write a TH1D of the eigenvalues
  TH1D *EigenVals_ratio = new TH1D("Eigen_Values_Ratio", "Eigen_Values_Ratio", npars, 0, npars);
  for (int i = 0; i < npars; ++i) {
    EigenVals_ratio->SetBinContent(i+1, (eigen->GetEigenValues())(i)/eigentot);
  }
  EigenVals_ratio->Write("Eigen_Values_ratio");

  // Write a TH1D of the eigenvalues
  double eigen_cum = 0;
  TH1D *EigenVals_cum = new TH1D("Eigen_Values_Cum", "Eigen_Values_Cum", npars, 0, npars);
  for (int i = 0; i < npars; ++i) {
    eigen_cum += (eigen->GetEigenValues())(i);
    EigenVals_cum->SetBinContent(i+1, eigen_cum/eigentot);
  }
  EigenVals_cum->Write("Eigen_Values_cum");

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

    temp_th2->Write(Form("Recreated_%i", transfer[i].GetNcols()));
    // Copy
    temp_th2->Divide(OrigTH2);
    temp_th2->SetMinimum(0.5);
    temp_th2->SetMaximum(1.5);
    temp_th2_copy->SetMinimum(-0.1);
    temp_th2_copy->SetMaximum(0.1);
    temp_th2->Write(Form("Recreated_%i_RelNom", transfer[i].GetNcols()));
    // Make percentage difference
    temp_th2_copy->Write(Form("Recreated_%i_RelNomPercent", transfer[i].GetNcols()));

    // Make a temp 1D object
    TH1D *parvals = new TH1D(Form("Parvals_%i", transfer[i].GetNcols()), Form("Parvals_%i;Parameter number;Parameter value", transfer[i].GetNcols()), Throws[i]->GetXaxis()->GetNbins(), 0, Throws[i]->GetXaxis()->GetNbins());
    TH1D *parerr = new TH1D(Form("Parerr_%i", transfer[i].GetNcols()), Form("Parerr_%i;Parameter number;Parameter uncertainty", transfer[i].GetNcols()), Throws[i]->GetXaxis()->GetNbins(), 0, Throws[i]->GetXaxis()->GetNbins());

    // Now make the 1D central value and error band
    for (int j = 0; j < npars; ++j) {
      // Project each of the 2D throws vs parameter onto the 1D axis, getting mean and rms
      TH1D *projecty = (TH1D*)Throws[i]->ProjectionY("_py", j+1, j+1);
      parvals->SetBinContent(j+1, projecty->GetMean());
      parerr->SetBinContent(j+1, projecty->GetRMS());
    }

    // Loop over bins and make relative the actual value
    TH1D *parerr_rel = (TH1D*)parerr->Clone(Form("%s_rel", parerr->GetName()));
    parerr_rel->GetYaxis()->SetTitle("Relative difference to actual uncertainty");
    for (int j = 0; j < npars; ++j) {
      parerr_rel->SetBinContent(j+1, (parerr_rel->GetBinContent(j+1)-sqrt((*origmat)(j,j)))/sqrt((*origmat)(j,j)));
    }

    parvals->Write();
    parerr->Write();
    parerr_rel->Write();
    // And write the throws
    Throws[i]->Write();
    delete parvals;
    delete parerr;
    delete parerr_rel;
  }

  // Close the file
  std::cout << "Wrote output to " << output->GetName() << std::endl;
  output->Close();
}

void makedecomp(std::string filename) {

  // The input file
  TFile *file = new TFile(filename.c_str());

  // The collection of matrices
  //TMatrixDSym *mat = (TMatrixDSym*)file->Get("Covariance_Matrix")->Clone();
  TMatrixDSym *mat = (TMatrixDSym*)file->Get("nddet_cov")->Clone();
  file->Close();

  // Create the large object
  MatrixCollect matobj(mat);

  // Construct all the new matrices
  matobj.ConstructNew();

  // Test the new matrices
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
