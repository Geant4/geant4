{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("testem2.root");
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist1 = (TH1D*)f.Get("1");
  hist1->Draw("HIST");
  
  TH1D* hist4 = (TH1D*)f.Get("4");
  hist4->Draw("HIST");
  
  TH1D* hist5 = (TH1D*)f.Get("5");
  hist5->Draw("HIST");  
}  
