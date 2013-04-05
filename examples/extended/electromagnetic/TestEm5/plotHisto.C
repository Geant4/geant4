{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("kulchi.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist12 = (TH1D*)f.Get("12");
  hist12->Draw("HIST");
  
  TH1D* hist13 = (TH1D*)f.Get("13");
  hist13->Draw("HIST");
}  
