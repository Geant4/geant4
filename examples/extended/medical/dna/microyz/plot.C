// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
//
// Author: Sebastien Incerti, CNRS, France
// Date: 25 Feb. 2015
// The Geant4-DNA collaboration
// *********************************************************************
{
gROOT->Reset();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
gStyle->SetOptStat(00000);

//***************************************
//*************************************** 
// MAKE YOUR SELECTIONS
// for histograms
//***************************************
//***************************************

Int_t linB=100;    // linear histo: nb of bins in x - 1000 is best for integration
Double_t ymin=1;   // minimum x-axis value
Double_t ymax=300; // maximum x-axis value

//***************************************
//***************************************

system ("rm -rf yz.root");
system ("hadd yz.root yz_*.root");

c1 = new TCanvas ("c1","",60,60,800,800);
Int_t mycolor;

TFile f("yz.root"); 
mycolor=4;

TNtuple* ntuple;
ntuple = (TNtuple*)f.Get("yz"); 

Double_t radius,eventID,nofHits,nbEdep,y,z,Einc;
ntuple->SetBranchAddress("radius",&radius);
ntuple->SetBranchAddress("eventID",&eventID);
ntuple->SetBranchAddress("nbHits",&nofHits);
ntuple->SetBranchAddress("nbScoredHits",&nbEdep);
ntuple->SetBranchAddress("y",&y);
ntuple->SetBranchAddress("z",&z);
ntuple->SetBranchAddress("Einc",&Einc);

//plot f(y)

c1->cd(1);

TH1F *hfyw = new TH1F ("hfyw","hfyw",linB,0,ymax);

Int_t nentries = (Int_t)ntuple->GetEntries();
Double_t population=0;
Double_t yLocalMin=1e100;
Double_t yLocalMax=0;

Double_t yF_anal=0;
Double_t yD_anal=0;

for (Int_t i=0; i<nentries; i++) 
{
 ntuple->GetEntry(i);

 hfyw->Fill(y,nofHits/nbEdep);
 if (yLocalMin>y) yLocalMin=y;
 if (yLocalMax<y) yLocalMax=y;
 population=population+nofHits/nbEdep;
 yF_anal = yF_anal + (nofHits/nbEdep)*y;
 yD_anal = yD_anal + (nofHits/nbEdep)*y*y;
}

cout << "**** Results ****" << endl;
cout << endl;
cout << "---> yF =" << yF_anal/population << " keV/um" << endl;
cout << "---> yD =" << (yD_anal/population)/(yF_anal/population) << " keV/um" << endl;
cout << endl;
cout << "---> Limits: " << endl;
cout << "     * min value of y = " << yLocalMin << " keV/um" << endl;
cout << "     * max value of y = " << yLocalMax << " keV/um" << endl;

if ( (yLocalMax>ymax) || (yLocalMin<ymin) ) 
{
  cout << "WARNING: please check your histogram limits ! " << endl;
}

gPad->SetLogy();
hfyw->Scale (1./(population*hfyw->GetBinWidth(1)));
hfyw->SetTitle("f(y) (um/keV)");
hfyw->GetXaxis()->SetTitle("y (keV/um)");
hfyw->SetFillColor(2);
hfyw->SetLineColor(2);
hfyw->Draw("HIST");
}

