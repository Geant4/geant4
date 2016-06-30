{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f("allproc.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");

  c1->Divide(3,1);

  c1->cd(1);
  TH1D* hist1 = (TH1D*)f.Get("h1");
  //  c1->SetLogy(1);
  hist1->Draw("HIST");

  TH1D* hist11 = (TH1D*)f.Get("h11");
  hist11->SetLineColor(2);
  hist11->Draw("HIST SAME");

  c1->cd(2);
  // c1->SetLogy(1);
  TH1D* hist2 = (TH1D*)f.Get("h2");
  hist2->Draw("HIST");
  
  TH1D* hist12 = (TH1D*)f.Get("h12");
  hist12->SetLineColor(2);
  hist12->Draw("HIST SAME");

  c1->cd(3);
  // c1->SetLogy(1);
  TH1D* hist3 = (TH1D*)f.Get("h3");
  hist3->Draw("HIST");
    
  TH1D* hist13 = (TH1D*)f.Get("h13");
  hist13->SetLineColor(2);
  hist13->Draw("HIST SAME");          

  c1->Print("allproc.gif");
  c1->Close();
}  
