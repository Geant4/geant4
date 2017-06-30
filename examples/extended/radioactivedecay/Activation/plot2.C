{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  c1->SetLogy(1);
  c1->cd();
  c1->Update();
  
  TFile fa("Li7.root");  
  TFile fb("Li6.root");  

  TH1D* ha1 = (TH1D*)fa.Get("6");  
  ha1->SetStats(kFALSE);
  ha1->SetLineColor(kBlue);        
  ha1->Draw("HIST");
  TH1D* hb1 = (TH1D*)fb.Get("6");  
  hb1->SetLineColor(kRed);   
  hb1->Draw("HIST SAME");
   
  // Print the histograms legend
  TLegend *le = new TLegend(0.6,0.6,0.8,0.8);
  le->AddEntry(ha1,"Li7","l");  
  le->AddEntry(hb1,"Li6","l");
  le->Draw();      
}  
