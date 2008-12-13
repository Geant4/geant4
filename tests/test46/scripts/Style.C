{ 
gROOT->Reset();
gROOT->SetStyle("Plain"); 
gStyle->SetOptStat(0);

Int_t nplot = 4;
Int_t nener = 3;
Int_t ndir  = 3;

Int_t iplot = 0;
Int_t iener = 0;
Int_t idir = 0;

TString axtit[3] = {"E (GeV)", "","log10(p/GeV)"};
TString tit[nplot] = {"ECAL energy deposit", "HCAL visible energy","HCAL energy deposit","Total energy deposit"};

TString hed[ndir] = {"QBBC","QBBC+ch.ex","QBBC+ch.ex+xsec"};
TString dir[ndir] = {"QBBC","QBBCF","QBBCG"};
TString fil[nener]= {"pi-8gev","pi-30gev","pi-300gev"};
TString part[nener]= {"#pi^{-} 8 GeV","#pi^{-} 30 GeV","#pi^{-} 300 GeV"};
Double_t ener[nener]= {8.,30.,300.};
Double_t x1[nplot] = {1.0,0.02,1.0,2.0}; 
//Double_t x2[nxs] = {5,6,5,6,6,6,6}; 
Double_t y1[nplot] = {0.32,0.1,0.16,0.3}; 
//Double_t y2[nxs] = {20,10,20,5,1,1,1}; 
TLegend* leg[nplot];
TH1F*    hh[nplot];
TString  gtit, fin;

TFile   ff[ndir];
Int_t   col[ndir]  = {2, 3, 4};

}
