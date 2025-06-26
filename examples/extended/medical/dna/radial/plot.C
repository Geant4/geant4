// -------------------------------------------------------------------
// -------------------------------------------------------------------
//
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

  TCanvas* c1 = new TCanvas ("c1","",20,20,800,800);
  c1->Divide(1,1);

  TFile* f = new TFile("radial.root");

  TNtuple* ntuple;
  ntuple = (TNtuple*)f->Get("radial");
  bool rowWise = true;
  TBranch* eventBranch = ntuple->FindBranch("row_wise_branch");
  if ( ! eventBranch ) rowWise = false;

  c1->cd(1);
  gStyle->SetOptStat(000000);
  gPad->SetLogy();

  ntuple->SetLineWidth(5);
  ntuple->Draw("dose:radius","","L");

  TH1* hist = (TH1*)gPad->GetPrimitive("htemp");
  hist->SetTitle("Absorbed dose VS radius");
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetTitleSize(0.03);
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->SetTitleOffset(1.8);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitle("Radius (nm)");
  hist->GetYaxis()->SetTitle("Dose (Gy)");
  gPad->SetTicks(1, 1);
  gPad->Modified();
}

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address) {
  TLeaf* leaf = ntuple->FindLeaf(name);
  if ( ! leaf ) {
    std::cerr << "Error in <SetLeafAddress>: unknown leaf --> " << name << std::endl;
    return;
  }
  leaf->SetAddress(address);
}
