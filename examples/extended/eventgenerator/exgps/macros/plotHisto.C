// Simple ROOT macro to plot histograms from test02.
// Only active histograms are processed.

{
  gROOT->Reset();

  // Draw histos filled by Geant4 simulation
  //
  TFile f = TFile("test02.root");
  auto c1 = new TCanvas("test02", "test02", 200, 10, 1200, 700);
  c1->Divide(4,2);

  auto dir = (TDirectory*)f.Get("histo");

  auto hist1 = (TH1D*)dir->Get("h1.1");
  c1->cd(1);
  hist1->Draw("HIST");

  auto hist3 = (TH1D*)dir->Get("h1.3");
  c1->cd(2);
  hist3->Draw("HIST");

  auto hist4 = (TH1D*)dir->Get("h1.4");
  c1->cd(3);
  hist4->Draw("HIST");

  auto hist5 = (TH2D*)dir->Get("h2.1");
  c1->cd(5);
  hist5->Draw("HIST");

  auto hist8 = (TH2D*)dir->Get("h2.4");
  c1->cd(6);
  hist8->Draw("HIST");
}
