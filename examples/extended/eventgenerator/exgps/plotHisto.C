{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f = TFile("test02.root");  
  TCanvas* c1 = new TCanvas("c1", "  ");
  
  TDirectory* dir = f.Get("histo");
    
  TH1D* hist1 = (TH1D*)dir->Get("h1.1");      
  hist1->Draw("HIST");
/*  
  TH1D* hist2 = (TH1D*)dir->Get("h1.2");      
  hist2->Draw("HIST");
*/  
  TH1D* hist3 = (TH1D*)dir->Get("h1.3");    
  hist3->Draw("HIST");
  
  TH1D* hist4 = (TH1D*)dir->Get("h1.4");      
  hist4->Draw("HIST");
  
  TH2D* hist5 = (TH2D*)dir->Get("h2.1");      
  hist5->Draw("HIST");
/*  
  TH2D* hist6 = (TH2D*)dir->Get("h2.2");      
  hist6->Draw("HIST");
    
  TH2D* hist7 = (TH2D*)dir->Get("h2.3");      
  hist7->Draw("HIST");
*/    
  TH2D* hist8 = (TH2D*)dir->Get("h2.4");      
  hist8->Draw("HIST");                
}  
