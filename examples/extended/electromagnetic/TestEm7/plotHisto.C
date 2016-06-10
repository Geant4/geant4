{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f("testem7.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist1 = (TH1D*)f.Get("h1");
  hist1->Draw("HIST");
  
  TH1D* hist2 = (TH1D*)f.Get("h2");
  hist2->Draw("HIST");
  
  TH1D* hist3 = (TH1D*)f.Get("h3");
  hist3->Draw("HIST");
}  
