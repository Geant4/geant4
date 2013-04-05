{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("pixe.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");
  c1->SetLogy(1);
  c1->cd();
  c1->Update();
     
  TH1D* hist1 = (TH1D*)f.Get("3");
  hist1->Draw("HIST");
     
  TH1D* hist2 = (TH1D*)f.Get("5");
  hist2->Draw("HIST");
    
  TH1D* hist3 = (TH1D*)f.Get("20");
  hist3->Draw("HIST");
  
  TH1D* hist4 = (TH1D*)f.Get("40");
  hist4->Draw("HIST");      
}  
