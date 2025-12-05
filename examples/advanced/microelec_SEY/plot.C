// *********************************************************************
// To execute this macro under ROOT, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// *********************************************************************
{
gROOT->Reset();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
Double_t scale;
	
c1 = new TCanvas ("c1","",20,20,1000,500);
c1->Divide(2,1);

//system ("rm -rf microelectronics.root");
//system ("hadd microelectronics.root microelectronics_*.root");

TFile f("microelectronics.root"); 

TNtuple* ntuple;
ntuple = (TNtuple*)f.Get("microelectronics"); 
     
c1->cd(1);
  gStyle->SetOptStat(000000);
  
  // All
  ntuple->Draw("flagProcess","","B");
  ntuple->SetFillColor(2);
  
  // Elastic
  ntuple->Draw("flagProcess","flagProcess==11","same");
  ntuple->SetFillColor(3); 
  
  // Ionisation
  ntuple->Draw("flagProcess","flagProcess==12||flagProcess==14||flagProcess==15||flagProcess==16||flagProcess==17","same");
  ntuple->SetFillColor(4);
  
  gPad->SetLogy();

c1->cd(2);

  // Electrons
  ntuple->SetMarkerColor(2);
  ntuple->Draw("x:y:z/1000","flagParticle==1");

  // Protons
  ntuple->SetMarkerColor(4);
  ntuple->SetMarkerSize(4);
  ntuple->Draw("x:y:z/1000","flagParticle==2","same");

  // Ions
  ntuple->SetMarkerColor(3);
  ntuple->SetMarkerSize(3);
  ntuple->Draw("x:y:z/1000","flagParticle==3","same");
  
}
