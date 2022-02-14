// -------------------------------------------------------------------
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// *********************************************************************

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address);

void plot()
{
<<<<<<< HEAD
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
=======
  gROOT->Reset();
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");
  	
  TCanvas* c1 = new TCanvas ("c1","",20,20,1500,500);
  c1->Divide(3,1);
  
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
       
  // canvas tab 1
  c1->cd(1);
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
  
  // canvas tab 2
  c1->cd(2);

  ntuple->SetMarkerColor(2);
  ntuple->Draw("x:y:z/1000","flagParticle==1");

  //ntuple->SetMarkerColor(4);
  //ntuple->SetMarkerSize(4);
  //ntuple->Draw("x:y:z/1000","flagParticle==4 || flagParticle==5 || flagParticle==6","same");

<<<<<<< HEAD
end:
=======
  // canvas tab 3
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

  if ( ! rowWise ) {
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
  }
  else {
    SetLeafAddress(ntuple, "flagParticle",&flagParticle);
    SetLeafAddress(ntuple, "flagProcess",&flagProcess);
    SetLeafAddress(ntuple, "x",&x);
    SetLeafAddress(ntuple, "y",&y);
    SetLeafAddress(ntuple, "z",&z);
    SetLeafAddress(ntuple, "totalEnergyDeposit",&totalEnergyDeposit);
    SetLeafAddress(ntuple, "stepLength",&stepLength);
    SetLeafAddress(ntuple, "kineticEnergyDifference",&kineticEnergyDifference);
    SetLeafAddress(ntuple, "kineticEnergy",&kineticEnergy);
    SetLeafAddress(ntuple, "cosTheta",&angle);
    SetLeafAddress(ntuple, "eventID",&eventID);
    SetLeafAddress(ntuple, "trackID",&trackID);
    SetLeafAddress(ntuple, "parentID",&parentID);
    SetLeafAddress(ntuple, "stepID",&stepID);
  }

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
}

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address) {
  TLeaf* leaf = ntuple->FindLeaf(name);
  if ( ! leaf ) {
    std::cerr << "Error in <SetLeafAddress>: unknown leaf --> " << name << std::endl;
    return;
  }
  leaf->SetAddress(address);
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}
