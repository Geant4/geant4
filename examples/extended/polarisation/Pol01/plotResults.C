
void plotResults(const char* fname = "polar.root") {
  // Open a ROOT file and read histograms
  TFile myfile(fname,"READ");
  TH1D* histo[12];

  for (Int_t j=0; j<12;j++)
    histo[j] = (TH1D *) myfile.Get(Form("h%d",j+1));

  // Produce a plot
  TCanvas *c1 = new TCanvas("c1","The Pol01 example",1200,900);
  gStyle->SetOptStat(0);
  c1->Divide(4,3);

  for (Int_t j=0; j<12;j++) {
    c1->cd(j+1);
    histo[j]->ls();
    histo[j]->DrawCopy();
  }
  c1->Update();
}
