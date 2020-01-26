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

int main(int argc, char** argv) {
  /*
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
      //mat(i,j) = corr(i,j);
    }
  }
  */

  if (argc != 2) {
    std::cerr << "Please give filename" << std::endl;
    return -1;
  }
  std::string filename = std::string(argv[1]);

  TFile *file = new TFile(filename.c_str());
  TMatrixDSym *mat = (TMatrixDSym*)file->Get("nddet_cov")->Clone();
  if (mat == NULL) {
    std::cerr << "Couldn't find nddet_cov in " << filename << std::endl;
    return -1;
  }
  int ndim = mat->GetNcols();

  //std::cout << "Original cov matrix:" << std::endl;
  //mat.Print();

  // Make the eigenvectors
  TMatrixDSymEigen eigen(*mat);
  TVectorD evals = eigen.GetEigenValues();
  TMatrixD evecs = eigen.GetEigenVectors();

  //std::cout << "Eigenvalues:" << std::endl;
  //evals.Print();
  //std::cout << "Eigenvectors:" << std::endl;
  //evecs.Print();

  // Start the PCA
  int nvals = ndim*1/2;
  TMatrixD transfer(evecs.GetSub(0, evecs.GetNcols()-1, 0, nvals-1));
  //std::cout << "Transfer matrix:" << std::endl;
  //transfer.Print();

  TMatrixD transferT(transfer);
  transferT.T();
  //std::cout << "Transfer matrix transpose:" << std::endl;
  //transferT.Print();

  // Make our parameter values
  TVectorD *pars = (TVectorD*)file->Get("det_weights");
  //pars(0) = 1.07;
  //pars(1) = 0.96;
  //pars(2) = 0.96;

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
  TDecompChol chol(*mat);
  chol.Decompose();
  TMatrixD cholmat = chol.GetU();
  cholmat.T();

  // Make throws
  const int nthrows = 1000;
  TH2D *throws = new TH2D("throws", "throws;parameter number; parameter value", ndim, 0, ndim, 100, 0, 2);
  TRandom3 *rand = new TRandom3(5);

  TH2D *throws_eigen = new TH2D("throws_eigen", "throws_eigen;parameter number; parameter value", ndim, 0, ndim, 100, 0, 2);

  // The random throws
  TVectorD randoms(ndim);
  // Now let's also make some throws in the eigen value basis and transform out into the actual parameter space
  TVectorD randoms_eigen(nvals);

  for (int i = 0; i < nthrows; ++i) {

    if (i % (nthrows/20) == 0) std::cout << "On throw " << i << "/" << nthrows << " (" << double(i)/nthrows*100 << "%)" << std::endl;
    for (int j = 0; j < ndim; ++j) {
      randoms(j) = rand->Gaus(0,1);
      if (j < nvals) randoms_eigen(j) = rand->Gaus(0,1);
    }
    /*
    std::cout << "Randoms in eigen basis:" << std::endl;
    randoms_eigen.Print();
    std::cout << "Transfer matrix: " << std::endl;
    transfer.Print();
    std::cout << "Transposed transfer matrix: " << std::endl;
    transferT.Print();
    */
    TVectorD rands_of_eigen = (transfer*randoms_eigen);
    //std::cout << "randoms_eigen*transfer:" << std::endl;
    //rands_eigen.Print();

    // X = AZ
    TVectorD rands = (cholmat)*randoms;
    TVectorD rands_chol_eigen = (cholmat)*rands_of_eigen;

    // X = E sqrt(eigen) Z
    TVectorD rands_eigen = evecs*fullsqrt*randoms;
    //(reduced*evecs).Print();
    //full.Print();
    //evecs.Print();
    //rands_eigen.Print();
    for (int j = 0; j < ndim; ++j) {
      throws->Fill(j, (*pars)(j)+rands(j));
      throws_eigen->Fill(j, (*pars)(j)+rands_chol_eigen(j));

      //throws_eigen->Fill(j, pars(j)+rands_eigen(j));
      //throws->Fill(j, rands(j));
      //throws_eigen->Fill(j, rands_eigen(j));
    }
  }

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  TString canvname = filename+Form("%i_pdf", nvals);
  canv->Print(canvname+".pdf[");
  canv->cd();
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);
  throws->SetTitle(Form("N_{throws}=%i", nthrows));
  throws_eigen->SetTitle(Form("N_{throws}=%i", nthrows));
  throws->Scale(1./nthrows);
  throws_eigen->Scale(1./nthrows);
  throws->Draw("colz");
  canv->Print(canvname+".pdf");
  canv->Clear();
  throws_eigen->Draw("colz");
  canv->Print(canvname+".pdf");
  TH1D *throw1d = new TH1D("throws1d", "throws1d; parameter number; parameter value", ndim, 0, ndim);
  TH1D *throw1d_eigen = new TH1D("throws1d_eigen", "throws1d_eigen; parameter number; parameter value", ndim, 0, ndim);

  for (int i = 0; i < ndim; ++i) {
    TH1D *projy = throws->ProjectionY("_py", i+1, i+1);
    TH1D *projy2 = throws_eigen->ProjectionY("_py", i+1, i+1);
    canv->Clear();
    projy->SetLineColor(kBlue);
    projy2->SetLineColor(kRed);
    double maximum = projy->GetMaximum() > projy2->GetMaximum() ? projy->GetMaximum() : projy2->GetMaximum();
    projy->Draw();
    projy2->Draw("same");
    projy->GetYaxis()->SetRangeUser(0, maximum*1.2);
    canv->Print(canvname+".pdf");
    //throw1d->SetBinContent(i+1, projy->GetMean());
    throw1d->SetBinContent(i+1, projy->GetRMS()-sqrt((*mat)(i,i)));
    //throw1d_eigen->SetBinContent(i+1, projy2->GetMean());
    throw1d_eigen->SetBinContent(i+1, projy2->GetRMS()-sqrt((*mat)(i,i)));
    delete projy;
    delete projy2;
  }

  throw1d->SetTitle("Residual");
  throw1d->GetYaxis()->SetTitle("Residual (estimate-truth)");
  canv->Clear();
  double maximum = throw1d->GetMaximum() > fabs(throw1d->GetMinimum()) ? throw1d->GetMaximum() : fabs(throw1d->GetMinimum());
  maximum = throw1d_eigen->GetMaximum() > maximum ? throw1d_eigen->GetMaximum() : maximum;
  maximum = fabs(throw1d_eigen->GetMinimum()) > maximum ? fabs(throw1d_eigen->GetMinimum()) : maximum;
  maximum *= 1.2;
  throw1d->SetMinimum(-1*maximum);
  throw1d->SetMaximum(maximum);
  throw1d->Draw();
  throw1d->SetLineColor(kBlue);
  throw1d_eigen->Draw("same");
  //truth->Draw("same");
  //throw1d_eigen->SetLineStyle(kDashed);
  throw1d_eigen->SetLineColor(kRed);
  canv->Print(canvname+".pdf");
  canv->Print(canvname+".pdf]");


  TFile *out = new TFile(Form("%s_%i_anal.root", filename.c_str(), nvals), "recreate");
  throws->Write();
  throws_eigen->Write();
  throw1d->Write();
  throw1d_eigen->Write();
  out->Close();

  return 0;
}
