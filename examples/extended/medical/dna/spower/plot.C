// *********************************************************************
// To execute this macro under ROOT after your simulation ended,
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// *********************************************************************

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address);

void plot()
{
  gROOT->Reset();
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");
  TCanvas* c1 = new TCanvas ("c1","",20,20,500,500);

  c1->Divide(1,1);

  TFile* f = new TFile("spower.root");

  TNtuple* ntuple;
  ntuple = (TNtuple*)f->Get("track");
  bool rowWise = true;
  TBranch* eventBranch = ntuple->FindBranch("row_wise_branch");
  if ( ! eventBranch ) rowWise = false;

  c1->cd(1);
  gStyle->SetOptStat(000000);
  //gPad->SetLogx();
  gPad->SetLogy();

  ntuple->SetFillStyle(1001);
  ntuple->SetFillColor(2);

  TH1F *htmp = new TH1F("htmp", "Secondary electron energy distribution per incident primary", 100, 0, 1000);

  ntuple->Draw("kineticEnergy>>htmp","kineticEnergy<100000&&flagParticle==1","B");

  htmp->Scale (1./ntuple->GetEntries("trackID==1"));
  htmp->SetXTitle("Energy (eV)");
  htmp->Draw();
}

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address) {
  TLeaf* leaf = ntuple->FindLeaf(name);
  if ( ! leaf ) {
    std::cerr << "Error in <SetLeafAddress>: unknown leaf --> " << name << std::endl;
    return;
  }
  leaf->SetAddress(address);
}
