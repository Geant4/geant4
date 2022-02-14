{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  c1->SetLogx(1);
  c1->SetLogy(1);
  c1->SetGridy(1);
  c1->cd();
  c1->Update();
  
  TFile fa("pixe.root");  

  TH1D* ha1 = (TH1D*)fa.Get("11");  
  ha1->SetLineColor(kBlue);        
  ha1->Draw("HIST");

}
