{ 
gROOT->Reset();
gROOT->SetStyle("Plain"); 
gStyle->SetOptStat(0);

Int_t ntarg = 5;
Int_t npart = 4;
Int_t nxs   = 6;

Int_t ipart = 0;
Int_t itarg = 0;
Int_t ixs   = 0;

TString axtit[3] = {"log10(E/MeV)", "Cross section (bn)","log10(p/GeV)"};

TString hed[nxs] = {"Elastic","Elastic","Inelastic","Inelastic","Capture","Fission"};
TString fil0[nxs]= {"elp","el","inlp","inl","cap","fis"};
Double_t x1[nxs] = {-3,-1,-3,0,-1,-1}; 
Double_t x2[nxs] = {5,6,5,6,6,6}; 
Double_t y1[nxs] = {0.001,0.001,0.001,0.001,0.1,0.1}; 
Double_t y2[nxs] = {20,10,20,5,10,10}; 
TLegend* leg[nxs];
TH1F*    hh[nxs];
TString  gtit, fin;

TFile   ff[ntarg];
TString targ[ntarg] = {"H","C","Al","Fe","Pb"};
TString filt[ntarg] = {"h","c","al","fe","pb"};
Int_t   col[ntarg]  = {2, 3, 4, 5, 6};

TString filp[npart] = {"p_","n_","pi+","pi-"};
TString part[npart] = {"proton","neutron","#pi^{+}","#pi^{-}"};

}
