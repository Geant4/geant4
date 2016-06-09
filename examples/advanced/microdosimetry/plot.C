// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// This macro needs the track.txt file
// *********************************************************************
{
gROOT->Reset();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
Double_t scale;
	
c1 = new TCanvas ("c1","",20,20,1000,500);
c1.Divide(2,1);

TFile f("microdosimetry.root"); 

TNtuple* ntuple;
ntuple = (TNtuple*)f->Get("ntuple"); 
     
c1.cd(1);
  gStyle->SetOptStat(000000);
  
  // All
  ntuple->Draw("flagProcess","","B");
  ntuple->SetFillColor(2);
  
  // Excitation
 
  ntuple->Draw("flagProcess","flagProcess==12||flagProcess==15||flagProcess==17||flagProcess==20||flagProcess==23||flagProcess==26||flagProcess==30","Bsame");
  ntuple->SetFillColor(3);
  
  // Elastic
  ntuple->Draw("flagProcess","flagProcess==11","Bsame");
  ntuple->SetFillColor(4);
  
  // Ionisation
  ntuple->Draw("flagProcess","flagProcess==13||flagProcess==18||flagProcess==21||flagProcess==24||flagProcess==27||flagProcess==31||flagProcess==33||flagProcess==34","Bsame");
  ntuple->SetFillColor(5);
  
  // Charge decrease
  ntuple->Draw("flagProcess","flagProcess==19||flagProcess==25||flagProcess==28","Bsame");
  ntuple->SetFillColor(6);
  
  // Charge increase
  ntuple->Draw("flagProcess","flagProcess==22||flagProcess==29||flagProcess==32","Bsame");
  
  gPad->SetLogy();

c1.cd(2);

  // Electrons
  ntuple->SetMarkerColor(2);
  ntuple->Draw("x:y:z/1000","flagParticle==1","");

  // Protons
  ntuple->SetMarkerColor(4);
  ntuple->SetMarkerSize(4);
  ntuple->Draw("x:y:z/1000","flagParticle==2","same");

  //Hydrogen
  ntuple->SetMarkerColor(3);
  ntuple->SetMarkerSize(3);
  ntuple->Draw("x:y:z/1000","flagParticle==3","same");
  
}
