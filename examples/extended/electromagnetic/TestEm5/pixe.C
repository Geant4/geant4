// ROOT macro file for plotting example B4 histograms 
// 
// Can be run from ROOT session:
// root[0] .x plotHisto.C

{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  
  // Draw histos filled by Geant4 simulation 
  //   

  // Open file filled by Geant4 simulation 
  TFile f("pixe.root");

  // Create a canvas and divide it into 2x2 pads
  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
  c1->Divide(2,2);
  
  // Draw pixe in the pad 1
  c1->cd(1);
  gPad->SetLogy(1);
  gPad->SetGridy(1);
  TH1D* hist1 = (TH1D*)f.Get("h55");
  hist1->SetLineColor(kBlue); 
  hist1->Draw("HIST");
  
  // Draw deexcitation in the pad 2
  c1->cd(2);
  gPad->SetLogy(1);
  gPad->SetGridy(0);  
  TH1D* hist2 = (TH1D*)f.Get("h51");
  hist2->SetLineColor(kBlue);   
  hist2->Draw("HIST");
  
  // Draw total in the pad 3
  c1->cd(3);
  gPad->SetLogy(1);
  gPad->SetGridy(1);
  TH1D* hist3 = (TH1D*)f.Get("h3");
  hist3->SetLineColor(kRed); 
  hist3->Draw("HIST");
}  
