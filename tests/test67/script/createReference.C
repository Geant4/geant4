#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TObjString.h"

int createReference()
{
  TFile* f = new TFile("reference_g4.root","RECREATE");
  /*
    KEY: TGraphErrors     PeakEff;1
    KEY: TGraphErrors     FullEff;1
    KEY: TObjString       Version;1       Collectable string class
    KEY: TObjString       List;1  Collectable string class
  */
  Int_t npoints = 6;
  TGraphErrors* gr1 = new TGraphErrors(npoints);
  gr1->SetName("PeakEff");
  TGraphErrors* gr2 = new TGraphErrors(npoints);
  gr2->SetName("FullEff");
  
  //
  gr1->SetPoint(0,45.,0.131);
  gr1->SetPointError(0,0.,0.001);
  gr1->SetPoint(1,60.,0.747);
  gr1->SetPointError(1,0.,0.002);
  gr1->SetPoint(2,120.,2.527);
  gr1->SetPointError(2,0.,0.016);
  gr1->SetPoint(3,200.,2.443);
  gr1->SetPointError(3,0.,0.003);
  gr1->SetPoint(4,500.,1.396);
  gr1->SetPointError(4,0.,0.004);
  gr1->SetPoint(5,2000.,0.576);
  gr1->SetPointError(5,0.,0.004);
  //
  gr2->SetPoint(0,45.,0.205);
  gr2->SetPointError(0,0.,0.001);
  gr2->SetPoint(1,60.,1.349);
  gr2->SetPointError(1,0.,0.004);
  gr2->SetPoint(2,120.,6.139);
  gr2->SetPointError(2,0.,0.018);
  gr2->SetPoint(3,200.,7.540);
  gr2->SetPointError(3,0.,0.009);
  gr2->SetPoint(4,500.,7.264);
  gr2->SetPointError(4,0.,0.012);
  gr2->SetPoint(5,2000.,5.448);
  gr2->SetPointError(5,0.,0.011);

  TObjString* stri = new TObjString("reference_g4");
  f->cd();
  gr1->Write(gr1->GetName());
  gr2->Write(gr2->GetName());
  stri->Write("Version");
  f->Close();

  return 0;
}
