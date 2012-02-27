{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("gammaSpectrum.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist1 = (TH1D*)f.Get("3");
  hist1->Draw("HIST");
  
  TH1D* hist2 = (TH1D*)f.Get("5");
  hist2->Draw("HIST");  
}  
