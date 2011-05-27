#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TObjString.h"
#include "TMath.h"
#include <iostream>
#include "TF1.h"

using namespace std;

int compare(TString file1,TString file2="reference_g4.root")
{
  gROOT->SetStyle("Plain");
  TFile* f1 = new TFile(file1);
  TFile* f2 = new TFile(file2);
  /*
    KEY: TGraphErrors     PeakEff;1
    KEY: TGraphErrors     FullEff;1
    KEY: TObjString       Version;1       Collectable string class
    KEY: TObjString       List;1  Collectable string class
  */
  TGraphErrors* gr1 = (TGraphErrors*) f1->Get("PeakEff");
  TGraphErrors* gr2 = (TGraphErrors*) f2->Get("PeakEff");

  Int_t npoints = gr1->GetN();
  if (gr2->GetN() != gr1->GetN())
    {
      std::cout << "The graphs have a different number of points!" << std::endl;
      return 1;
    }
  TGraphErrors* grdiff = new TGraphErrors(npoints);
  for (size_t i=0;i<npoints;i++)
    {
      Double_t x1,y1,x2,y2;
      gr1->GetPoint(i,x1,y1);
      gr2->GetPoint(i,x2,y2);
      Double_t e1 = gr1->GetErrorY(i);
      Double_t e2 = gr2->GetErrorY(i);
      if (TMath::Abs(x1-x2)>0.01*x1)
	{
	  cout << "The graphs have a different x grid!" << endl;
	  cout << x1 << " " << x2 << endl;
	  return 1;
	}
      grdiff->SetPoint(i,x1,y2-y1); 
      Double_t errdiff = TMath::Sqrt(e1*e1+e2*e2);
      grdiff->SetPointError(i,0.,errdiff);
    }
  grdiff->Draw("AZP");
  TF1* pol0 = new TF1("pol0","pol0(0)",0,2500);
  grdiff->Fit(pol0,"","",0,2500);
  cout << pol0->GetChisquare() << "/" << pol0->GetNDF() << " " << pol0->GetProb() << endl;
  //Force zero, to check the chi2
  //pol0->SetParameter(0,0.);
  //grdiff->Fit(pol0,"","",0,2500);
  //cout << pol0->GetChisquare() << "/" << pol0->GetNDF() << endl;

  return 0;
}
