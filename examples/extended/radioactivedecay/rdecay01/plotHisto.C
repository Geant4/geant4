{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("rdecay01.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  c1->SetLogy(0);
  c1->cd();
  c1->Update();
     
  TH1D* hist3 = (TH1D*)f.Get("3");
  hist3->Draw("HIST");    
}  
