// *********************************************************************
// To execute this macro under ROOT,
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// This macro needs the output ROOT file
// *********************************************************************

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address);

void plot()
{
  gROOT->Reset();
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");
  
  TCanvas* c1 = new TCanvas ("c1","",20,20,1000,500);
  
  TFile* f = new TFile("ICSD.root");
  f->ls();
  
  TH1D* hist1 = new TH1D("histo","ICSD", 30, 0.0, 30);

  Double_t ion = 0.0;

  TNtuple* ntuple1 = (TNtuple*)f->Get("ntuple_1");
  bool rowWise = true;
  TBranch* eventBranch = ntuple1->FindBranch("row_wise_branch");
  if ( ! eventBranch ) rowWise = false;
  std::cout <<  "rowWise: " << rowWise << std::endl; 

  if ( ! rowWise ) {
    ntuple1->SetBranchAddress("ionisations", &ion);
  } else {
    SetLeafAddress(ntuple1, "ionisations",&ion);
  }

  Int_t nentries = ntuple1->GetEntries();

  for (Int_t i=0; i<nentries; i++)
  {
      ntuple1->GetEntry(i);
      hist1->Fill(ion);
  }

  hist1->Draw();
  hist1->GetXaxis()->SetLabelSize(0.025);
  hist1->GetYaxis()->SetLabelSize(0.025);

  hist1->GetXaxis()->SetTitleSize(0.035);
  hist1->GetYaxis()->SetTitleSize(0.035);

//  hist1->GetXaxis()->SetTittleOffset(1.4);
//  hist1->GetYaxis()->SetTittleOffset(1.4);

  hist1->GetXaxis()->SetTitle("ionisation number");
  hist1->GetYaxis()->SetTitle("frequency");

  c1->SaveAs("ICSD.tiff");
}

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address) {
  TLeaf* leaf = ntuple->FindLeaf(name);
  if ( ! leaf ) {
    std::cerr << "Error in <SetLeafAddress>: unknown leaf --> " << name << std::endl;
    return;
  }
  leaf->SetAddress(address);
}
