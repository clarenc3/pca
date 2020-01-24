#include <iostream>
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDBase.h"
#include "TDecompChol.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"

int main() {
  // Make the original covariance matrix
  // Just a jokey one from TN315
  const int ndim = 3;
  TMatrixDSym corr(ndim);
  corr(0,0) = 1.00;
  corr(0,1) = -0.83;
  corr(0,2) = -0.01;
  corr(1,0) = -0.83;
  corr(1,1) = 1.00;
  corr(1,2) = -0.31;
  corr(2,0) = -0.01;
  corr(2,1) = -0.31;
  corr(2,2) = 1.00;

  // Now make the correlation matrix to a covariance matrix
  TVectorD uncert(ndim);
  uncert(0) = 0.15;
  uncert(1) = 0.15;
  uncert(2) = 0.40;

  TMatrixDSym mat(ndim);
  for (int i = 0; i < ndim; ++i) {
    for (int j = 0; j < ndim; ++j) {
      mat(i,j) = corr(i,j)*uncert(i)*uncert(j);
    }
  }

  //std::cout << "Original cov matrix:" << std::endl;
  //mat.Print();

  // Make the eigenvectors
  TMatrixDSymEigen eigen(mat);
  TVectorD evals = eigen.GetEigenValues();
  TMatrixD evecs = eigen.GetEigenVectors();

  //std::cout << "Eigenvalues:" << std::endl;
  //evals.Print();
  //std::cout << "Eigenvectors:" << std::endl;
  //evecs.Print();

  // Start the PCA
  const int nvals = 2;
  TMatrixD transfer(evecs.GetSub(0, evecs.GetNcols()-1, 0, nvals-1));
  //std::cout << "Transfer matrix:" << std::endl;
  //transfer.Print();

  TMatrixD transferT(transfer);
  transferT.T();
  //std::cout << "Transfer matrix transpose:" << std::endl;
  //transferT.Print();

  // Make our parameter values
  TVectorD pars(ndim);
  pars(0) = 1.07;
  pars(1) = 0.96;
  pars(2) = 0.96;

  // The original transfer
  //TMatrixD evecsT(evecs);
  //evecsT.T();
  //std::cout << "Full transfer:" << std::endl;
  //(evecsT*pars).Print();

  //std::cout << "Partial transfer:" << std::endl;
  //(transferT*pars).Print();

  TMatrixD fullsqrt(ndim, ndim);
  for (int i = 0; i < ndim; ++i) fullsqrt(i,i) = sqrt(evals(i));

  // Make the reduced eigen val matrix
  TMatrixD reduced(nvals, nvals);
  for (int i = 0; i < nvals; ++i) reduced(i,i) = evals(i);
  //std::cout << "Reduced eigen: " << std::endl;
  //reduced.Print();

  // New corr matrix
  //std::cout << "Recomposed matrix:" << std::endl;
  TMatrixD recomposed = (transfer*reduced*transferT);
  //recomposed.Print();

  //std::cout << "Original matrix:" << std::endl;
  //mat.Print();

  //std::cout << "Cholesky of original: " << std::endl;
  TDecompChol chol(mat);
  chol.Decompose();
  //chol.Print();


  // Make throws
  const int nthrows = 10000000;
  TH2D *throws = new TH2D("throws", "throws;parameter number; parameter value", ndim, 0, ndim, 1000, -1, 3);
  TRandom3 *rand = new TRandom3(5);

  TH2D *throws_eigen = new TH2D("throws_eigen", "throws_eigen;parameter number; parameter value", ndim, 0, ndim, 1000, -1, 3);

  // The random throws
  TVectorD randoms(ndim);
  for (int i = 0; i < nthrows; ++i) {
    for (int j = 0; j < ndim; ++j) {
      randoms(j) = rand->Gaus(0,1);
    }
    // X = AZ
    TVectorD rands = chol.GetU()*randoms;

    // X = E sqrt(eigen) Z
    TVectorD rands_eigen = evecs*fullsqrt*randoms;
    //(reduced*evecs).Print();
    //full.Print();
    //evecs.Print();
    //rands_eigen.Print();
    for (int j = 0; j < ndim; ++j) {
      throws->Fill(j, pars(j)+rands(j));
      throws_eigen->Fill(j, pars(j)+rands_eigen(j));
    }
  }

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->Print("asdf.pdf[");
  canv->cd();
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);
  throws->SetTitle(Form("N_{throws}=%i", nthrows));
  throws_eigen->SetTitle(Form("N_{throws}=%i", nthrows));
  throws->Scale(1./nthrows);
  throws_eigen->Scale(1./nthrows);
  throws->Draw("colz");
  canv->Print("asdf.pdf");
  canv->Clear();
  throws_eigen->Draw("colz");
  canv->Print("asdf.pdf");
  TH1D *throw1d = new TH1D("throws1d", "throws1d; parameter number; parameter value", ndim, 0, ndim);
  for (int i = 0; i < ndim; ++i) {
    TH1D *projy = throws->ProjectionY("_py", i+1, i+1);
    canv->Clear();
    projy->Draw();
    canv->Print("asdf.pdf");
    throw1d->SetBinContent(i+1, projy->GetMean());
    throw1d->SetBinError(i+1, projy->GetRMS());
  }
  canv->Clear();
  throw1d->Draw();
  canv->Print("asdf.pdf");
  canv->Print("asdf.pdf]");

  TFile *file = new TFile("asdf.root", "recreate");
  throws->Write();
  throw1d->Write();
  throws_eigen->Write();
  file->Close();


  return 0;
}
