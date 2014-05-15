{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("run1.root");
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist1 = (TH1D*)f.Get("h1");
  hist1->Draw("HIST");
  
  TH1D* hist2 = (TH1D*)f.Get("h2");
  hist2->Draw("HIST");
  
  TH1D* hist3 = (TH1D*)f.Get("h3");
  c1->SetLogy(1);
  c1->cd();
  c1->Update(); 
  hist3->Draw("HIST");
  
  TH1D* hist4 = (TH1D*)f.Get("h4");
  hist4->Draw("HIST");
  
  TH1D* hist5 = (TH1D*)f.Get("h5");
  hist5->Draw("HIST");
  
  TH1D* hist6 = (TH1D*)f.Get("h6");
  hist6->Draw("HIST");                  
}  
