{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("rdecay02.root");
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  c1->SetLogy(0);
  c1->cd();
  c1->Update();
  
  ///TH1D* hist1 = (TH1D*)f.Get("H11");
  ///hist1->SetLineColor(kRed); 
  ///hist1->Draw("HIST");
  
  ///TH1D* hist2 = (TH1D*)f.Get("H12");
  ///hist2->Draw("HIST");
     
  ///TH1D* hist3 = (TH1D*)f.Get("H13");
  ///hist3->SetLineColor(kRed);   
  ///hist3->Draw("HIST");
  
  ///TH1D* hist4 = (TH1D*)f.Get("H14");
  ///hist4->Draw("HIST");
  
  ///TH1D* hist5 = (TH1D*)f.Get("H15");
  ///hist5->Draw("HIST");
  
  TH1D* hist6 = (TH1D*)f.Get("H16");
  hist6->Draw("HIST");        
}  
