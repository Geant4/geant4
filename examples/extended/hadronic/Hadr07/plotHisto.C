{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f("bone.root");
  TCanvas* c1 = new TCanvas("c1", "  ");
  c1->SetLogy(0);
  c1->cd();
  c1->Update();
     
  TH1D* hist10 = (TH1D*)f.Get("10");
  hist10->Draw("HIST");

  TH1D* hist11 = (TH1D*)f.Get("11");
  hist11->Draw("HIST");
    
  TH1D* hist12 = (TH1D*)f.Get("12");
  hist12->Draw("HIST");
  
  TH1D* hist13 = (TH1D*)f.Get("13");
  hist13->Draw("HIST");

}  
