// -------------------------------------------------------------------
// $Id: plot.C,v 1.6 2010-10-08 10:01:35 sincerti Exp $
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

TNtuple* ntuple0;
ntuple0 = (TNtuple*)f->Get("ntuple0"); 
     
c1.cd(1);
  gStyle->SetOptStat(000000);
  
  // All
  ntuple0->Draw("flagProcess","");
  ntuple0->SetFillColor(2);
  
  // Excitation
 
  ntuple0->Draw("flagProcess","flagProcess==12||flagProcess==15||flagProcess==17||flagProcess==20||flagProcess==23||flagProcess==26||flagProcess==30","same");
  ntuple0->SetFillColor(3);
  
  // Elastic
  ntuple0->Draw("flagProcess","flagProcess==11","same");
  ntuple0->SetFillColor(4);
  
  // Ionisation
  ntuple0->Draw("flagProcess","flagProcess==13||flagProcess==18||flagProcess==21||flagProcess==24||flagProcess==27||flagProcess==31||flagProcess==33||flagProcess==34","same");
  ntuple0->SetFillColor(5);
  
  // Charge decrease
  ntuple0->Draw("flagProcess","flagProcess==19||flagProcess==25||flagProcess==28","same");
  ntuple0->SetFillColor(6);
  
  // Charge increase
  ntuple0->Draw("flagProcess","flagProcess==22||flagProcess==29||flagProcess==32","same");
  
  gPad->SetLogy();

c1.cd(2);

  // Electrons
  ntuple0->SetMarkerColor(2);
  ntuple0->Draw("x:y:z/1000","flagParticle==1");

  // Protons
  ntuple0->SetMarkerColor(4);
  ntuple0->SetMarkerSize(4);
  ntuple0->Draw("x:y:z/1000","flagParticle==2","same");

  //Hydrogen
  ntuple0->SetMarkerColor(3);
  ntuple0->SetMarkerSize(3);
  ntuple0->Draw("x:y:z/1000","flagParticle==3","same");
  
}
