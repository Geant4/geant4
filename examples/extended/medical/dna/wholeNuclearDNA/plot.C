// *********************************************************************
// To execute this macro under ROOT, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// This macro needs the output ROOT file
// *********************************************************************
{
gROOT->Reset();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
Double_t scale;
	
c1 = new TCanvas ("c1","",20,20,1000,500);

// The command below stands for:
// IF 1 THEN 2
system ("ls wholeNuclearDNA_* 2> /dev/null && hadd -f wholeNuclearDNA.root wholeNuclearDNA_*.root");

TFile f("wholeNuclearDNA.root"); 

TNtuple* ntuple;
ntuple = (TNtuple*)f.Get("ntuple"); 
     
c1->cd(1);
gStyle->SetOptStat(000000);
  
ntuple->SetMarkerColor(2);
ntuple->SetMarkerStyle(20);
ntuple->SetMarkerSize(1);
ntuple->Draw("x/1000:y/1000:z/1000","flagVolume==1","");

ntuple->SetMarkerColor(4);
ntuple->SetMarkerStyle(24);
ntuple->Draw("x/1000:y/1000:z/1000","flagVolume==2","same");
  
}
