{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f("run01.root");
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist1 = (TH1D*)f.Get("h1");
  hist1->Draw("HIST");
  
  TH1D* hist4 = (TH1D*)f.Get("h4");
  hist4->Draw("HIST");
  
  TH1D* prof4 = (TH1D*)f.Get("p4");
  prof4->Draw(" ");
    
  TH1D* hist5 = (TH1D*)f.Get("h5");
  hist5->Draw("HIST");
  
  TH1D* hist8 = (TH1D*)f.Get("h8");
  hist8->Draw("HIST");
  
  TH1D* prof8 = (TH1D*)f.Get("p8");
  prof8->Draw(" ");  
}  
