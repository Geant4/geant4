{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f("pixe.root");
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist3 = (TH1D*)f.Get("3");
  hist3->Draw("HIST");
  
  c1->SetLogy(1);
  c1->cd();
  c1->Update();
     
  TH1D* hist5 = (TH1D*)f.Get("5");
  hist5->Draw("HIST");
}  
