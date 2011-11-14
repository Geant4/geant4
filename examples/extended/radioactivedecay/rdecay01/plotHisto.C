{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("rdecay1.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  c1->SetLogy(1);
  c1->cd();
  c1->Update();
     
  TH1D* hist3 = (TH1D*)f.Get("6");
  hist3->Draw("HIST");    
}  
