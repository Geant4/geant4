{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f("essai.root");      
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  ///TH1D* hist1 = (TH1D*)f.Get("1");
  ///hist1->Draw("HIST");
  
  TH1D* hist2 = (TH1D*)f.Get("2");
  hist2->Draw("HIST");  
     
  TH1D* hist6 = (TH1D*)f.Get("6");
  hist6->Draw("HIST");
}  
