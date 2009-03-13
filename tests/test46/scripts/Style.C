{ 
gROOT->Reset();
gROOT->SetStyle("Plain"); 
gStyle->SetOptStat(0);

Int_t nplot = 4;
//Int_t nener = 17;
Int_t nener = 3;
Int_t ndir  = 7;

Int_t iplot = 0;
Int_t iener = 0;
Int_t idir = 0;

TString axtit[3] = {"E (GeV)", "","log10(p/GeV)"};
TString tit[nplot] = {"ECAL energy deposit", "HCAL visible energy","HCAL energy deposit","Total energy deposit"};

TString hed[3] = {"QBBC","QBBC+ch.ex","QBBC+ch.ex+xsec"};
TString dir[ndir] = {"QGSP_BERT_EMV","QGSP_BERT","FTFP_BERT","FTFP_BERT_EML","QBBC","QBBCF","QBBCG"};
TString fil[nener]= {"pi-9gev","pi-30gev","pi-300gev"};
TString part[nener]= {"#pi^{-} 9 GeV","#pi^{-} 30 GeV","#pi^{-} 300 GeV"};
//TString part[nener]= {"proton 9 GeV","proton 30 GeV","proton 300 GeV"};
Double_t ener[nener]= {9.,30.,300.};
//Double_t ener[nener]= {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,20.,30.,50.,100.,200.,300.};
Double_t x1[nplot] = {1.0,0.015,1.5,1.5}; 
//Double_t x2[nxs] = {5,6,5,6,6,6,6}; 
Double_t y1[nplot] = {0.12,0.06,0.04,0.20}; 
//Double_t y2[nxs] = {20,10,20,5,1,1,1}; 
TLegend* leg[nplot];
TH1F*    hh[nplot];
TString  gtit, fin;

TFile   ff[ndir];
Int_t   col[ndir]  = {2, 3, 4, 5, 6, 7, 1};
Int_t   iac[ndir]  = {0, 0, 0, 0, 0, 0, 0};

}
