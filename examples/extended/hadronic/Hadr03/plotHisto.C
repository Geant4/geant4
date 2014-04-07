{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  ////TFile f = TFile("elastic.root");
  TFile f = TFile("Li6a.root");    
  TCanvas* c1 = new TCanvas("c1", "  ");
  c1->SetLogy(1);
  c1->cd();
  c1->Update();   
  
  ///TH1D* hist1 = (TH1D*)f.Get("1");
  ///hist1->Draw("HIST");
  
  ///TH1D* hist2 = (TH1D*)f.Get("2");
  ///hist2->Draw("HIST");
  
  TH1D* hist3 = (TH1D*)f.Get("3");
  hist3->Draw("HIST");
  
  TH1D* hist4 = (TH1D*)f.Get("4");
  hist4->Draw("HIST");
  
  TH1D* hist5 = (TH1D*)f.Get("5");
  hist5->Draw("HIST");
      
  TH1D* hist11 = (TH1D*)f.Get("11");
  hist11->Draw("HIST");
  
  ///TH1D* hist12 = (TH1D*)f.Get("12");
  ///hist12->Draw("HIST");    
}  
