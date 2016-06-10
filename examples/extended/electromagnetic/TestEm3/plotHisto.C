{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  ///TFile f("testem3.root");
  TFile f("run01.root");    
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TH1D* hist1 = (TH1D*)f.Get("1");
  hist1->Draw("HIST");
  
  TH1D* hist2 = (TH1D*)f.Get("2");
  hist2->Draw("HIST");
  
  TH1D* hist3 = (TH1D*)f.Get("11");
  c1->SetLogy(1);
  c1->cd();
  c1->Update(); 
  hist3->Draw("HIST");
  
  TH1D* hist4 = (TH1D*)f.Get("12");
  hist4->Draw("HIST");
  
  TH1D* hist5 = (TH1D*)f.Get("21");
  hist5->Draw("HIST");
  
  TH1D* hist6 = (TH1D*)f.Get("22");
  hist6->Draw("HIST");
          
}  
