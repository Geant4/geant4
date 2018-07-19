{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f("NeutronSource.root");      
  TCanvas* c1 = new TCanvas("c1", "  ");

  ///TH1D* hist4 = (TH1D*)f.Get("4");
  ///hist4->Draw("HIST");

  TH1D* hist6 = (TH1D*)f.Get("6");
  hist6->Draw("HIST");
}
