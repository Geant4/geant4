// -------------------------------------------------------------------
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT after your simulation ended,
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// *********************************************************************

void plot()
{
  gROOT->Reset();
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");

  TCanvas* c1 = new TCanvas ("c1","",20,20,1000,500);
  c1->Divide(2,1);

  // Uncomment if merging should be done
  //system ("rm -rf dna.root");
  //system ("hadd dna.root dna_*.root");

  TFile* f = new TFile("dna.root");

  TNtuple* ntuple;
  ntuple = (TNtuple*)f->Get("dna");
  bool rowWise = true;
  TBranch* eventBranch = ntuple->FindBranch("row_wise_branch");
  if ( ! eventBranch ) rowWise = false;
  // std::cout <<  "rowWise: " << rowWise << std::endl;

  // Canvas tab 1
  c1->cd(1);
  gStyle->SetOptStat(000000);

  // All
  ntuple->SetFillStyle(1001);
  ntuple->SetFillColor(2);
  ntuple->Draw("flagProcess","","B");

  // Excitation
  ntuple->SetFillStyle(1001);
  ntuple->SetFillColor(3);
  ntuple->Draw("flagProcess","flagProcess==12||flagProcess==15||flagProcess==22||flagProcess==32||flagProcess==42||flagProcess==52||flagProcess==62","Bsame");

  // Elastic
  ntuple->SetFillStyle(1001);
  ntuple->SetFillColor(4);
  ntuple->Draw("flagProcess","flagProcess==11||flagProcess==21||flagProcess==31||flagProcess==41||flagProcess==51||flagProcess==61||flagProcess==110||flagProcess==210||flagProcess==410||flagProcess==510||flagProcess==710||flagProcess==120||flagProcess==220||flagProcess==420||flagProcess==520||flagProcess==720","Bsame");

  // Ionisation
  ntuple->SetFillStyle(1001);
  ntuple->SetFillColor(5);
  ntuple->Draw("flagProcess","flagProcess==13||flagProcess==23||flagProcess==33||flagProcess==43||flagProcess==53||flagProcess==63||flagProcess==73||flagProcess==130||flagProcess==230||flagProcess==430||flagProcess==530||flagProcess==730","Bsame");

  // Charge decrease
  //ntuple->SetFillStyle(1001);
  //ntuple->SetFillColor(6);
  //ntuple->Draw("flagProcess","flagProcess==24||flagProcess==44||flagProcess==54","Bsame");

  // Charge increase
  //ntuple->SetFillStyle(1001);
  //ntuple->SetFillColor(7);
  //ntuple->Draw("flagProcess","flagProcess==35||flagProcess==55||flagProcess==65","Bsame");

  gPad->SetLogy();

  // Canvas tab 2
  c1->cd(2);

  // Electrons
  ntuple->SetMarkerColor(2);
  ntuple->SetMarkerStyle(20);
  ntuple->SetMarkerSize(.2);
  ntuple->Draw("x:y:z","flagParticle==1");

  // Protons, hydrogen
  ntuple->SetMarkerColor(4);
  ntuple->Draw("x:y:z","flagParticle==2 || flagParticle==3 ","same");
}
