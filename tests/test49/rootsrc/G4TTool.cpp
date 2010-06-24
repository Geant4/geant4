//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TTool.h"


ClassImp(G4TTool)


using namespace std;
using namespace ROOT;
using namespace TMath;


//______________________________________________________________________________
void G4TTool::Initialize()
{
  fHxmin = 3;
  fHxmax = fPublication->fHeader.fTypeValue + 3;
  fHymin = 1.0001e-8;
  fHymax = 2.0001e-8;
  fHfbin =  1;
  fHnbin =  50;
  fDangle =  9.0; // @@ It should not be a parameter! It should be taken from DB!
  fPi = 3.14159265;
  fDanrad = fPi * fDangle / 180;
  fHlxmin = std::log(fHxmin);
  fHlxmax = std::log(fHxmax);
  //PrepareHistograms(fHnbin, fHlxmin, fHlxmax); // It should be done, when xmax is defined
}

//______________________________________________________________________________
void G4TTool::PrepareHistograms(Int_t hnbin, Double_t hlxmin, Double_t hlxmax)
{
  // Check if deleted
  fHZL = (TH1F*)gDirectory->Get("hZL");
  fHDT = (TH1F*)gDirectory->Get("hDT");
  if(fHZL) delete fHZL;
  if(fHDT) delete fHDT;
  // create zero level and dT
  fHZL = new TH1F("hZL", "Zero Level", hnbin, hlxmin, hlxmax);
  fHDT = new TH1F("hDT", "dT in MeV", hnbin, hlxmin, hlxmax);
  // adjust zero levels, create dT histogram
  Double_t hzero =  0.1e-12;
  Double_t halfbn = 0.5 * (hlxmax - hlxmin) / hnbin;
  cout<<"===>>G4TTool::PrepareHistograms: hlxmax="<<hlxmax<<"("<<std::exp(hlxmax)
       <<"), hlxmin="<<hlxmin<<", hnbin="<<hnbin<<", halfbn="<<halfbn<<endl;
  Double_t* vZeroLevels = new Double_t[(Int_t)(hnbin+1)]; // 0 is underflow
  Double_t* pZeroLevels = new Double_t[(Int_t)(hnbin+1)];
  Double_t* vErrors = new Double_t[(Int_t)(hnbin+1)];
  for(Int_t i = 0; i<= hnbin; ++i)
  {
    vZeroLevels[i] = hzero;
    vErrors[i] = vZeroLevels[i] * 0.8;
  }
  // set the content and error
  fHZL->SetContent(&vZeroLevels[0]);
  fHZL->SetError(&vErrors[0]);
  // get values of center of bins abscissa into vector
  TAxis* axis = fHDT->GetXaxis();
  for(Int_t i = 0; i<= hnbin; ++i) pZeroLevels[i] = axis->GetBinCenter(i);
  // debug print
  Double_t sum=fHxmin;
  //
  vZeroLevels[0] = 1.;
  for(Int_t i = 1; i<= hnbin; ++i)
  {
    vZeroLevels[i] = std::exp(pZeroLevels[i] + halfbn) - std::exp(pZeroLevels[i] - halfbn);
    // debug print
    sum+=vZeroLevels[i];
    cout<<"G4TTool::PrepareHistograms:fHDT["<<i<<"] = "<<pZeroLevels[i]<<", d = "
	  <<vZeroLevels[i]<<", S = "<<sum<<endl;
    //
  }
  fHDT->SetContent(&vZeroLevels[0]);
  // clear
  delete vZeroLevels;
  delete pZeroLevels;
  delete vErrors;
}

//__________________________________________________________________________________
void G4TTool::RenderHSolid(TH1F* hist, Int_t hf, Int_t hn, Double_t m, Color_t color,
                           Bool_t noDots )
{
  if(hist)
  {
    cout<<"G4TTool::RenderHSolid: m="<<m<<endl;
    Int_t hnbin = hn;
    Double_t* vWx = new Double_t[hnbin];
    Double_t* vWy = new Double_t[hnbin];
    Double_t* vWd = new Double_t[hnbin];
    Double_t* vWt = new Double_t[hnbin];
    Double_t* vWp = new Double_t[hnbin];
    Double_t* vWf = new Double_t[hnbin];
    // get values of center of bins abscissa into vector
    TAxis* axis = hist->GetXaxis();
    if(!axis) cout<<"Warning>G4TTool::RenderHSolid: No axis found!"<<endl;
    for(Int_t i = 0; i < hn; ++i) // 0 is still underflow !
    {
      vWx[i] = axis->GetBinCenter(i+1);         // log(T)
      vWy[i] = hist->GetBinContent(i+1);        // N
    }
    for(Int_t i =  0; i < hnbin; ++i)
    {
      vWt[i] = std::exp(vWx[i]);                // T
      vWp[i] = sqrt(vWt[i] * (m + m + vWt[i])); // p
      vWf[i] = vWy[i] / vWp[i];                 // N/p
      cout<<"G4TTool::RenderHSolid: i#"<<i<<", T="<<vWt[i]<<", N="<<vWy[i]<<", p="<<vWp[i]
          <<", N/p="<<vWf[i]<<endl;
    }
    TGraph* graph = new TGraph(hnbin, vWt, vWf); 
    graph->SetLineColor(color); // red
    graph->SetLineWidth(2);
    graph->Draw();
    if(!noDots)
    {
      TGraph* dots = new TGraph(hnbin,vWt, vWf);
      dots->SetLineStyle(4); // dot dot
      dots->SetLineWidth(2);
      dots->Draw();
    }
    // clean up
    delete[] vWx;
    delete[] vWy;
    delete[] vWd;
    delete[] vWp;
    delete[] vWt;
    delete[] vWf;
  }
  else cout<< "Error>G4TTool::RenderHSolid: no histogram found!"<<endl;
}
