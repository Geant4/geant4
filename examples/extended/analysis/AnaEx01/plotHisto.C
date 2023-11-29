{
  gROOT->Reset();

  // Draw histos filled by Geant4 simulation
  //
  TFile f = TFile("AnaEx01.root");
  TCanvas* c1 = new TCanvas("c1", "  ");

  TDirectory* dir = f.Get("histo");

  TH1D* hist1 = (TH1D*)dir->Get("1");
  hist1->Draw("HIST");

  TH1D* hist2 = (TH1D*)dir->Get("2");
  hist2->Draw("HIST");

  TH1D* hist3 = (TH1D*)dir->Get("3");
  hist3->Draw("HIST");

  TH1D* hist4 = (TH1D*)dir->Get("4");
  hist4->Draw("HIST");
}
