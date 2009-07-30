{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("testem5.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist1 = (TH1D*)f.Get("13");
  hist1->Draw("HIST");
  
  c1->SetLogy(1);
  c1->cd();
  c1->Update(); 
  hist1->Draw("HIST");    
}  
