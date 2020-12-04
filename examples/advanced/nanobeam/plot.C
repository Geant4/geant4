{
// Use this macro for the visualisation of coefficient computation

  gROOT->Reset();
  gROOT->SetStyle("Plain");

  gStyle->SetOptStat(1111);
  gStyle->SetPalette(1);

  c1 = new TCanvas ("c1","",200,10,800,800);
  c1->Divide(2,2);

  system ("rm -rf nanobeam.root");
  system ("hadd -O nanobeam.root nanobeam_*.root");
  
  TFile f("nanobeam.root"); 

  TNtuple* ntuple0;
  ntuple0 = (TNtuple*)f.Get("ntuple0"); 
  
  c1->cd(1);
  ntuple0->SetMarkerSize(.2);
  ntuple0->SetMarkerColor(2);
  ntuple0->SetMarkerStyle(20);
  ntuple0->Draw("xIn:yIn:zIn","");
 
  c1->cd(2); 
  TNtuple* ntuple2;
  ntuple2 = (TNtuple*)f.Get("ntuple2"); 
  gStyle->SetPalette(1);
  ntuple2->SetMarkerColor(2);
  ntuple2->SetMarkerStyle(20);
  gPad->SetLogz();
  ntuple2->Draw("yIn:xIn","","hcolz");
 
  c1->cd(3);
  gStyle->SetPalette(1);
  ntuple2->SetMarkerColor(4);
  ntuple2->SetMarkerStyle(20);
  gPad->SetLogz();
  ntuple2->Draw("thetaIn:xIn","","hcolz");

  c1->cd(4);
  gStyle->SetPalette(1);
  ntuple2->SetMarkerColor(4);
  ntuple2->SetMarkerStyle(20);
  gPad->SetLogz();
  ntuple2->Draw("phiIn:yIn","","hcolz");
   
jump2:
/*
  c1->cd(1);
  gStyle->SetPalette(1);
  TNtuple* ntuple1;
  ntuple1 = (TNtuple*)f->Get("ntuple1"); 
  ntuple1->SetMarkerSize(.2);
  ntuple1->SetMarkerColor(2);
  ntuple1->SetMarkerStyle(20);
  gPad->SetLogz();
  ntuple1->Draw("yIn:xIn","","hcolz");
  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.4);
  htemp->GetYaxis()->SetTitleOffset(1.4);
  htemp->GetXaxis()->SetTitle("X (mm)");
  htemp->GetYaxis()->SetTitle("Y (mm)");
  htemp->SetTitle("Grid shadow");
*/ 
}
