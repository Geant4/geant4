{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  c1->SetLogy(0);
  c1->cd();
  c1->Update();
  
  TFile fa("Co60.root");  

  TH1D* ha1 = (TH1D*)fa.Get("4");  
  ha1->SetStats(kFALSE);
  ha1->SetLineColor(kRed);        
  ha1->Draw("HIST");

}
