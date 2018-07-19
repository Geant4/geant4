{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f("Co60.root");      
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  ///TH1D* hist1 = (TH1D*)f.Get("1");
  ///hist1->Draw("HIST");
  
  TH1D* hist4 = (TH1D*)f.Get("4");
  hist4->Draw("HIST");  
     
  TH1D* hist14 = (TH1D*)f.Get("14");
  hist14->Draw("HIST");
  
  TH1D* hist24 = (TH1D*)f.Get("24");
  hist24->Draw("HIST");  
     
  TH1D* hist25 = (TH1D*)f.Get("25");
  hist25->Draw("HIST");
  
  TH1D* hist26 = (TH1D*)f.Get("26");
  hist26->Draw("HIST");  
}  
