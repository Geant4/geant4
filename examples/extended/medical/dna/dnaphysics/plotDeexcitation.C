// -------------------------------------------------------------------
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT after your simulation ended,
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// *********************************************************************

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address);

void plotDeexcitation()
{
  gROOT->Reset();
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");

  TCanvas* c1 = new TCanvas ("c1","",20,20,1000,500);
  c1->Divide(2,1);

  // Uncomment if merging should be done
  //system ("rm -rf dna.root");
  //system ("hadd dna.root dna_*.root");

  TFile* f = new TFile("dna.root");

  TNtuple* ntuple2;
  ntuple2 = (TNtuple*)f->Get("track");
  bool rowWise2 = true;
  TBranch* eventBranch2 = ntuple2->FindBranch("row_wise_branch");
  if ( ! eventBranch2 ) rowWise2 = false;

  // Auger electrons
  c1->cd(1);
  gStyle->SetOptStat(000000);

  ntuple2->SetFillStyle(1001);
  ntuple2->SetFillColor(2);
  ntuple2->Draw("kineticEnergy","kineticEnergy>450&&kineticEnergy<550&&flagParticle==1","B");

  gPad->SetLogy();

  TText *pt1 = new TText(510,500,"Auger electrons");
  pt1->Draw("SAME");

  // Fluorescence photons
  c1->cd(2);
  gStyle->SetOptStat(000000);

  ntuple2->SetFillStyle(1001);
  ntuple2->SetFillColor(4);
  ntuple2->Draw("kineticEnergy","flagParticle==0","B");

  gPad->SetLogy();

  TText *pt2 = new TText(523.1,10,"Fluo. photons");
  pt2->Draw("SAME");

}

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address) {
  TLeaf* leaf = ntuple->FindLeaf(name);
  if ( ! leaf ) {
    std::cerr << "Error in <SetLeafAddress>: unknown leaf --> " << name << std::endl;
    return;
  }
  leaf->SetAddress(address);
}
