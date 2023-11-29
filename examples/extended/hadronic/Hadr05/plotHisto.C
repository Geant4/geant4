{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  ///TFile f("hadr05.root");
  TFile f("Fe-Sci.root");    
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist1 = (TH1D*)f.Get("1");
  hist1->Draw("HIST");
  
  TH1D* hist2 = (TH1D*)f.Get("2");
  hist2->Draw("HIST");
  
  TH1D* hist11 = (TH1D*)f.Get("11");
  c1->SetLogy(0);
  c1->cd();
  c1->Update(); 
  hist11->Draw("HIST");
  
  TH1D* hist12 = (TH1D*)f.Get("12");
  hist12->Draw("HIST");
  
  TH1D* hist21 = (TH1D*)f.Get("21");
  hist21->Draw("HIST");
  
  TH1D* hist22 = (TH1D*)f.Get("22");
  hist22->Draw("HIST");
  
  TH1D* hist23 = (TH1D*)f.Get("23");
  hist23->Draw("HIST");  
  
  TH1D* hist24 = (TH1D*)f.Get("24");
  hist24->Draw("HIST");              
}  
