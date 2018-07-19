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
	
c1 = new TCanvas ("c1","",20,20,1500,500);
c1->Divide(3,1);

system ("rm -rf dna.root");
system ("hadd dna.root dna_*.root");

TFile f("dna.root"); 

TNtuple* ntuple;
ntuple = (TNtuple*)f.Get("dna"); 
     
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
  
c1->cd(2);

  ntuple->SetMarkerColor(2);

  ntuple->Draw("x:y:z","flagParticle==1");

  //ntuple->SetMarkerColor(4);
  //ntuple->SetMarkerSize(4);
  //ntuple->Draw("x:y:z/1000","flagParticle==4 || flagParticle==5 || flagParticle==6","same");

c1->cd(3);

  Double_t flagParticle;
  Double_t flagProcess;
  Double_t x;
  Double_t y;
  Double_t z;
  Double_t totalEnergyDeposit;
  Double_t stepLength;
  Double_t kineticEnergyDifference;
  Int_t eventID;
  Double_t kineticEnergy;
  Int_t stepID;
  Int_t trackID;
  Int_t parentID;
  Double_t angle;

  ntuple->SetBranchAddress("flagParticle",&flagParticle);
  ntuple->SetBranchAddress("flagProcess",&flagProcess);
  ntuple->SetBranchAddress("x",&x);
  ntuple->SetBranchAddress("y",&y);
  ntuple->SetBranchAddress("z",&z);
  ntuple->SetBranchAddress("totalEnergyDeposit",&totalEnergyDeposit);
  ntuple->SetBranchAddress("stepLength",&stepLength);
  ntuple->SetBranchAddress("kineticEnergyDifference",&kineticEnergyDifference);
  ntuple->SetBranchAddress("kineticEnergy",&kineticEnergy);
  ntuple->SetBranchAddress("cosTheta",&angle);
  ntuple->SetBranchAddress("eventID",&eventID);
  ntuple->SetBranchAddress("trackID",&trackID);
  ntuple->SetBranchAddress("parentID",&parentID);
  ntuple->SetBranchAddress("stepID",&stepID);

  TH1F* hsolvE = new TH1F ("hsolvE","solvE",100,0,2000);
  TH1F* helastE = new TH1F ("helastE","elastE",100,0,2000);
  TH1F* hexcitE = new TH1F ("hexcitE","excitE",100,0,2000);
  TH1F* hioniE = new TH1F ("hiioniE","ioniE",100,0,2000);
  TH1F* hattE = new TH1F ("hattE","attE",100,0,2000);
  TH1F* hvibE = new TH1F ("hvibE","vibE",100,0,2000);
 
  for (Int_t j=0;j<ntuple->GetEntries(); j++) 
  {

    ntuple->GetEntry(j);
    if (flagProcess==10) hsolvE->Fill(x);
    if (flagProcess==11) helastE->Fill(x);
    if (flagProcess==12) hexcitE->Fill(x);
    if (flagProcess==13) hioniE->Fill(x);
    if (flagProcess==14) hattE->Fill(x);
    if (flagProcess==15) hvibE->Fill(x);

  }

  helastE->GetXaxis()->SetTitle("x (nm)");
  helastE->SetLineColor(2);

  hexcitE->SetLineColor(3);
  hioniE->SetLineColor(4);
  hattE->SetLineColor(5);
  hvibE->SetLineColor(6);
  hsolvE->SetLineColor(7);

  gPad->SetLogy();

  helastE->Draw("");
  hexcitE->Draw("SAME");
  hioniE->Draw("SAME");
  hattE->Draw("SAME");
  hvibE->Draw("SAME");
  hsolvE->Draw("SAME");


end:
}
