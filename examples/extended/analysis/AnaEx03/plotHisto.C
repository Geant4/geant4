{
  gROOT->Reset();

  // Draw histos filled by Geant4 simulation
  //
  TFile f = TFile("e-.root");
  TCanvas* c1 = new TCanvas("c1", "  ");
  c1->Divide(2,2);

  TDirectory* dir = (TDirectory*)f.Get("histo");

  c1->cd(1);
  TH1D* hist1 = (TH1D*)dir->Get("EAbs");
  hist1->Draw("HIST");

  c1->cd(2);
  TH1D* hist2 = (TH1D*)dir->Get("EGap");
  hist2->Draw("HIST");

  c1->cd(3);
  TH1D* hist3 = (TH1D*)dir->Get("LAbs");
  hist3->Draw("HIST");

  c1->cd(4);
  TH1D* hist4 = (TH1D*)dir->Get("LGap");
  hist4->Draw("HIST");
}
