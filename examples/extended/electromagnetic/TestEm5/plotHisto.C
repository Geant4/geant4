{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("athr.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist1 = (TH1D*)f.Get("1");
  hist1->Draw("HIST");
}  
