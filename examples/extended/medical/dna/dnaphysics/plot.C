// -------------------------------------------------------------------
// $Id: plot.C 70323 2013-05-29 07:57:44Z gcosmo $
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// *********************************************************************
{
gROOT->Reset();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
	
c1 = new TCanvas ("c1","",20,20,1000,500);
c1->Divide(2,1);

system ("rm -rf dna.root");
system ("hadd dna.root dna_*.root");

TFile f("dna.root"); 

TNtuple* ntuple;
ntuple = (TNtuple*)f.Get("dna"); 
     
c1->cd(1);
  gStyle->SetOptStat(000000);
  
  // All
  ntuple->Draw("flagProcess","","B");
  ntuple->SetFillStyle(1001);
  ntuple->SetFillColor(2);

  // Excitation
  ntuple->Draw("flagProcess","flagProcess==12||flagProcess==15||flagProcess==16||flagProcess==19||flagProcess==22||flagProcess==25||flagProcess==29","Bsame");
  ntuple->SetFillStyle(1001);
  ntuple->SetFillColor(3);

  // Elastic
  ntuple->Draw("flagProcess","flagProcess==11","Bsame");
  ntuple->SetFillStyle(1001);
  ntuple->SetFillColor(4);
  
  // Ionisation
  ntuple->Draw("flagProcess","flagProcess==13||flagProcess==17||flagProcess==20||flagProcess==23||flagProcess==26||flagProcess==30","Bsame");
  ntuple->SetFillStyle(1001);
  ntuple->SetFillColor(5);
  
  // Charge decrease
  //ntuple->Draw("flagProcess","flagProcess==18||flagProcess==24||flagProcess==27","Bsame");
  //ntuple->SetFillStyle(1001);
  //ntuple->SetFillColor(6);

  // Charge increase
  //ntuple->Draw("flagProcess","flagProcess==21||flagProcess==28||flagProcess==31","Bsame");
  //ntuple->SetFillStyle(1001);
  //ntuple->SetFillColor(7);
  
  gPad->SetLogy();
  
c1->cd(2);

  ntuple->SetMarkerColor(2);
  ntuple->Draw("x:y:z/1000","flagParticle==1");

  //ntuple->SetMarkerColor(4);
  //ntuple->SetMarkerSize(4);
  //ntuple->Draw("x:y:z/1000","flagParticle==4 || flagParticle==5 || flagParticle==6","same");

end:
}
