{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  ////TFile f = TFile("Hadr04.root");
  TFile f = TFile("Water_nothermal.root");      
  TCanvas* c1 = new TCanvas("c1", "  ");
  c1->SetLogy(0);
  c1->cd();
  c1->Update();   
  
  TH1D* hist1 = (TH1D*)f.Get("1");
  hist1->Draw("HIST");
  
  TH1D* hist2 = (TH1D*)f.Get("2");
  hist2->Draw("HIST");
  
  TH1D* hist3 = (TH1D*)f.Get("3");
  hist3->Draw("HIST");
  
  TH1D* hist4 = (TH1D*)f.Get("4");
  hist4->Draw("HIST");
  
  TH1D* hist5 = (TH1D*)f.Get("5");
  hist5->Draw("HIST");
  
  TH1D* hist6 = (TH1D*)f.Get("6");
  hist6->Draw("HIST");
  
  TH1D* hist7 = (TH1D*)f.Get("7");
  hist7->Draw("HIST");      
}  
