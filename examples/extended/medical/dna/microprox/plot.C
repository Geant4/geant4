// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
//
// Author: Sebastien Incerti
// Date: March 2nd, 2019
// The Geant4-DNA collaboration
// *********************************************************************

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address);

void plot()
{
gROOT->Reset();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
gStyle->SetOptStat(00000);

//********************
Int_t nbRadius = 101;
//********************

TCanvas *c1 = new TCanvas ("c1","",60,60,800,800);
Int_t mycolor;

TFile f("t.root"); 
mycolor=4;

TNtuple* ntuple;
ntuple = (TNtuple*)f.Get("t"); 

bool rowWise = true;
TBranch* eventBranch = ntuple->FindBranch("row_wise_branch");
if ( ! eventBranch ) rowWise = false;
// std::cout <<  "rowWise: " << rowWise << std::endl

Double_t radius1,nofHits,nbEdep,edep,radius2,Einc;
Int_t noRadius;

if ( ! rowWise ) {
  ntuple->SetBranchAddress("radius1",&radius1);
  ntuple->SetBranchAddress("noRadius",&noRadius);
  ntuple->SetBranchAddress("nbHits",&nofHits);
  ntuple->SetBranchAddress("nbScoredHits",&nbEdep);
  ntuple->SetBranchAddress("edep",&edep);
  ntuple->SetBranchAddress("radius2",&radius2);
  ntuple->SetBranchAddress("Einc",&Einc);
  }
else {
  SetLeafAddress(ntuple, "radius1",&radius1);
  SetLeafAddress(ntuple, "noRadius",&noRadius);
  SetLeafAddress(ntuple, "nbHits",&nofHits);
  SetLeafAddress(ntuple, "nbScoredHits",&nbEdep);
  SetLeafAddress(ntuple, "edep",&edep);
  SetLeafAddress(ntuple, "radius2",&radius2);
  SetLeafAddress(ntuple, "Einc",&Einc);        
}

Int_t nentries = (Int_t)ntuple->GetEntries();

//

Double_t t[1000]; // 1000 is the max number of radius values
Double_t population[1000];
Double_t myRad[1000]; 

for (Int_t i=0; i<1000; i++) 
{
 t[i]=0; 
 population[i]=0;
 myRad[i]=0; 
}

Int_t event = 0;

for (Int_t i=0; i<nentries; i++) 
{
  ntuple->GetEntry(i);
  t[noRadius] = t[noRadius] + edep;
  population[noRadius]=population[noRadius]+1;
  myRad[noRadius] = radius1;  
}

// Mean

for (Int_t j=1; j<nbRadius; j++) 
{
  t[j] = t[j]/population[j];
  t[j] = t[j]/(myRad[j+1]-myRad[j]);
  //cout << j << " " << myRad[j] << " " << myRad[j+1] 
  // << " " << t[j] << " " << population[j] << endl;
}

//

 c1->cd(1);

 TGraph* gr1 =new TGraph(nbRadius,myRad,t);
 gr1->SetMarkerColor(2);
 gr1->SetMarkerStyle(20);
 gr1->SetMarkerSize(1);
 gr1->SetLineColor(2);
 gr1->SetTitle("");
 gr1->GetXaxis()->SetLimits(0.1,100);
 gr1->GetYaxis()->SetLimits(0.,22);
 gr1->GetXaxis()->SetLabelSize(0.025);
 gr1->GetYaxis()->SetLabelSize(0.025);
 gr1->GetXaxis()->SetTitleSize(0.035);
 gr1->GetYaxis()->SetTitleSize(0.035);
 gr1->GetXaxis()->SetTitleOffset(1.4);
 gr1->GetYaxis()->SetTitleOffset(1.4);
 gr1->GetXaxis()->SetTitle("r (nm)");
 gr1->GetYaxis()->SetTitle("t (eV/nm)");
 gr1->Draw("");
 gPad->SetLogx();

}

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address) {
  TLeaf* leaf = ntuple->FindLeaf(name);
  if ( ! leaf ) {
    std::cerr << "Error in <SetLeafAddress>: unknown leaf --> " << name << std::endl;
    return;
  }
  leaf->SetAddress(address);
}

