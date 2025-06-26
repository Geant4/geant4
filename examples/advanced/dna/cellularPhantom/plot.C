// -------------------------------------------------------------------
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT,
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// Written by S. Incerti, 25/01/2025
// *********************************************************************
{
gROOT->Reset();
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0000);
gStyle->SetPalette(1);

auto c1 = new TCanvas ("c1","",20,20,1200,900);
c1->Divide(4,3);

//------------------------------
// Original phantom file view
//------------------------------

FILE * fp = fopen("phantoms/phantom.dat","r");

Double_t X, Y, Z, mat, tmp;
char unit[100];
Double_t voxelSizeX, voxelSizeY, voxelSizeZ;
Long_t numberVoxTot, numberVoxRed, numberVoxGreen, numberVoxBlue;

TNtuple *ntuplePhantom = new TNtuple("PHANTOM","ntuple","X:Y:Z:mat");

Long_t nlines=0;
Long_t ncols=0;

while (1)
   {
      if ( nlines == 0 ) ncols = fscanf(fp,"%ld %ld %ld %ld",&numberVoxTot,&numberVoxRed,&numberVoxGreen,&numberVoxBlue);
      if ( nlines == 1 ) ncols = fscanf(fp,"%lf %lf %lf %s",&tmp,&tmp,&tmp,unit);
      if ( nlines == 2 ) ncols = fscanf(fp,"%lf %lf %lf %s",&voxelSizeX,&voxelSizeY,&voxelSizeZ, unit);
      if ( nlines >= 3 ) ncols = fscanf(fp,"%lf %lf %lf %lf", &X, &Y, &Z, &mat);
      //cout << X << " " << Y << " " << Z << " " << mat << endl;
      if (ncols < 0) break;
      ntuplePhantom->Fill(X,Y,Z,mat);
      nlines++;
   }
fclose(fp);

c1->cd(1);

ntuplePhantom->SetMarkerColor(1);
ntuplePhantom->Draw("Y:X");
// RED
ntuplePhantom->SetMarkerColor(2);
ntuplePhantom->Draw("Y:X","mat==1","same");
// GREEN
ntuplePhantom->SetMarkerColor(3);
ntuplePhantom->Draw("Y:X","mat==2","same");
// BLUE
ntuplePhantom->SetMarkerColor(4);
ntuplePhantom->Draw("Y:X","mat==3","same");
//
TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
htemp->GetXaxis()->SetTitle("X (microns)");
htemp->GetYaxis()->SetTitle("Y (mirons)");
htemp->GetXaxis()->SetLabelSize(0.025);
htemp->GetYaxis()->SetLabelSize(0.025);
htemp->GetXaxis()->SetTitleSize(0.035);
htemp->GetYaxis()->SetTitleSize(0.035);
htemp->GetXaxis()->SetTitleOffset(1.4);
htemp->GetYaxis()->SetTitleOffset(1.4);
htemp->SetTitle("RGB phantom YX view");

c1->cd(5);

ntuplePhantom->SetMarkerColor(1);
ntuplePhantom->Draw("Y:Z");
// RED
ntuplePhantom->SetMarkerColor(2);
ntuplePhantom->Draw("Y:Z","mat==1","same");
// GREEN
ntuplePhantom->SetMarkerColor(3);
ntuplePhantom->Draw("Y:Z","mat==2","same");
// BLUE
ntuplePhantom->SetMarkerColor(4);
ntuplePhantom->Draw("Y:Z","mat==3","same");
//
TH2F *htempBis = (TH2F*)gPad->GetPrimitive("htemp");
htempBis->GetXaxis()->SetTitle("Z (microns)");
htempBis->GetYaxis()->SetTitle("Y (mirons)");
htempBis->GetXaxis()->SetLabelSize(0.025);
htempBis->GetYaxis()->SetLabelSize(0.025);
htempBis->GetXaxis()->SetTitleSize(0.035);
htempBis->GetYaxis()->SetTitleSize(0.035);
htempBis->GetXaxis()->SetTitleOffset(1.4);
htempBis->GetYaxis()->SetTitleOffset(1.4);
htempBis->SetTitle("RGB phantom YZ view");

c1->cd(9);

ntuplePhantom->SetMarkerColor(1);
ntuplePhantom->Draw("X:Z");
// RED
ntuplePhantom->SetMarkerColor(2);
ntuplePhantom->Draw("X:Z","mat==1","same");
// GREEN
ntuplePhantom->SetMarkerColor(3);
ntuplePhantom->Draw("X:Z","mat==2","same");
// BLUE
ntuplePhantom->SetMarkerColor(4);
ntuplePhantom->Draw("X:Z","mat==3","same");
//
TH2F *htempTer = (TH2F*)gPad->GetPrimitive("htemp");
htempTer->GetXaxis()->SetTitle("Z (microns)");
htempTer->GetYaxis()->SetTitle("X (mirons)");
htempTer->GetXaxis()->SetLabelSize(0.025);
htempTer->GetYaxis()->SetLabelSize(0.025);
htempTer->GetXaxis()->SetTitleSize(0.035);
htempTer->GetYaxis()->SetTitleSize(0.035);
htempTer->GetXaxis()->SetTitleOffset(1.4);
htempTer->GetYaxis()->SetTitleOffset(1.4);
htempTer->SetTitle("RGB phantom XZ view");

//------------------
// Read ROOT file
//------------------

TFile *f = new TFile ("phantom.root");

TNtuple* ntuple1;
TNtuple* ntuple2;
TNtuple* ntuple3;

ntuple1 = (TNtuple*)f->Get("ntuple1");
ntuple2 = (TNtuple*)f->Get("ntuple2");
ntuple3 = (TNtuple*)f->Get("ntuple3");

Double_t x, y, z, energy, dose;
Int_t voxelID;

ntuple1->SetBranchAddress("x",&x);
ntuple1->SetBranchAddress("y",&y);
ntuple1->SetBranchAddress("z",&z);
ntuple1->SetBranchAddress("energy",&energy);
ntuple1->SetBranchAddress("dose",&dose);
ntuple1->SetBranchAddress("voxelID",&voxelID);

ntuple2->SetBranchAddress("x",&x);
ntuple2->SetBranchAddress("y",&y);
ntuple2->SetBranchAddress("z",&z);
ntuple2->SetBranchAddress("energy",&energy);
ntuple2->SetBranchAddress("dose",&dose);
ntuple2->SetBranchAddress("voxelID",&voxelID);

ntuple3->SetBranchAddress("x",&x);
ntuple3->SetBranchAddress("y",&y);
ntuple3->SetBranchAddress("z",&z);
ntuple3->SetBranchAddress("energy",&energy);
ntuple3->SetBranchAddress("dose",&dose);
ntuple3->SetBranchAddress("voxelID",&voxelID);

//---------------------------------
// Absorbed energy distributions
//---------------------------------

c1->cd(2);
gPad->SetLogy();
ntuple1->Draw("energy","energy>0");
TH1F *htemp2 = (TH1F*)gPad->GetPrimitive("htemp");
htemp2->GetXaxis()->SetTitle("Energy (keV)");
htemp2->GetXaxis()->SetLabelSize(0.025);
htemp2->GetXaxis()->SetTitleSize(0.035);
htemp2->GetXaxis()->SetTitleOffset(1.4);
htemp2->SetTitle("RED voxel energy");
htemp2->SetFillStyle(1001);
htemp2->SetFillColor(2);

c1->cd(6);
gPad->SetLogy();
ntuple2->Draw("energy","energy>0");
TH1F *htemp3 = (TH1F*)gPad->GetPrimitive("htemp");
htemp3->GetXaxis()->SetTitle("Energy (keV)");
htemp3->GetXaxis()->SetLabelSize(0.025);
htemp3->GetXaxis()->SetTitleSize(0.035);
htemp3->GetXaxis()->SetTitleOffset(1.4);
htemp3->SetTitle("GREEN voxel energy");
htemp3->SetFillStyle(1001);
htemp3->SetFillColor(3);

c1->cd(10);
gPad->SetLogy();
ntuple3->Draw("energy","energy>0");
TH1F *htemp4 = (TH1F*)gPad->GetPrimitive("htemp");
htemp4->GetXaxis()->SetTitle("Energy (keV)");
htemp4->GetXaxis()->SetLabelSize(0.025);
htemp4->GetXaxis()->SetTitleSize(0.035);
htemp4->GetXaxis()->SetTitleOffset(1.4);
htemp4->SetTitle("BLUE voxel energy");
htemp4->SetFillStyle(1001);
htemp4->SetFillColor(4);

//------------------------------
// Map of energy distribution
//------------------------------

c1->cd(3);
TH2F *histNrjRed = new TH2F("histNrjRed","histNrjRed",100,0,800,100,0,800);
ntuple1->Draw("y:x>>histNrjRed","energy","contz");
gPad->SetLogz();
histNrjRed->Draw("contz");
histNrjRed->GetXaxis()->SetTitle("X (microns)");
histNrjRed->GetYaxis()->SetTitle("Y (mirons)");
histNrjRed->GetZaxis()->SetTitle("Energy (keV)");
histNrjRed->GetXaxis()->SetLabelSize(0.025);
histNrjRed->GetYaxis()->SetLabelSize(0.025);
histNrjRed->GetZaxis()->SetLabelSize(0.025);
histNrjRed->GetXaxis()->SetTitleSize(0.035);
histNrjRed->GetYaxis()->SetTitleSize(0.035);
histNrjRed->GetZaxis()->SetTitleSize(0.035);
histNrjRed->GetXaxis()->SetTitleOffset(1.4);
histNrjRed->GetYaxis()->SetTitleOffset(1.4);
histNrjRed->GetZaxis()->SetTitleOffset(.6);
histNrjRed->SetTitle("Energy map for RED voxels");

c1->cd(7);
TH2F *histNrjGreen = new TH2F("histNrjGreen","histNrjGreen",100,0,800,100,0,800);
ntuple2->Draw("y:x>>histNrjGreen","energy","contz");
gPad->SetLogz();
histNrjGreen->Draw("contz");
histNrjGreen->GetXaxis()->SetTitle("X (microns)");
histNrjGreen->GetYaxis()->SetTitle("Y (mirons)");
histNrjGreen->GetZaxis()->SetTitle("Energy (keV)");
histNrjGreen->GetXaxis()->SetLabelSize(0.025);
histNrjGreen->GetYaxis()->SetLabelSize(0.025);
histNrjGreen->GetZaxis()->SetLabelSize(0.025);
histNrjGreen->GetXaxis()->SetTitleSize(0.035);
histNrjGreen->GetYaxis()->SetTitleSize(0.035);
histNrjGreen->GetZaxis()->SetTitleSize(0.035);
histNrjGreen->GetXaxis()->SetTitleOffset(1.4);
histNrjGreen->GetYaxis()->SetTitleOffset(1.4);
histNrjGreen->GetZaxis()->SetTitleOffset(.6);
histNrjGreen->SetTitle("Energy map for GREEN voxels");

c1->cd(11);
TH2F *histNrjBlue = new TH2F("histNrjBlue","histNrjBlue",100,0,800,100,0,800);
ntuple3->Draw("y:x>>histNrjBlue","energy","contz");
gPad->SetLogz();
histNrjBlue->Draw("contz");
histNrjBlue->GetXaxis()->SetTitle("X (microns)");
histNrjBlue->GetYaxis()->SetTitle("Y (mirons)");
histNrjBlue->GetZaxis()->SetTitle("Energy (keV)");
histNrjBlue->GetXaxis()->SetLabelSize(0.025);
histNrjBlue->GetYaxis()->SetLabelSize(0.025);
histNrjBlue->GetZaxis()->SetLabelSize(0.025);
histNrjBlue->GetXaxis()->SetTitleSize(0.035);
histNrjBlue->GetYaxis()->SetTitleSize(0.035);
histNrjBlue->GetZaxis()->SetTitleSize(0.035);
histNrjBlue->GetXaxis()->SetTitleOffset(1.4);
histNrjBlue->GetYaxis()->SetTitleOffset(1.4);
histNrjBlue->GetZaxis()->SetTitleOffset(.6);
histNrjBlue->SetTitle("Energy map for BLUE voxels");

//----------------------------
// Map of dose distribution
//----------------------------

c1->cd(4);
TH2F *histDoseRed = new TH2F("histDoseRed","histDoseRed",100,0,800,100,0,800);
// WARNING : dose scaling to mGy
ntuple1->Draw("y:x>>histDoseRed","dose/1000","contz");
//gPad->SetLogz();
histDoseRed->Draw("contz");
histDoseRed->GetXaxis()->SetTitle("X (microns)");
histDoseRed->GetYaxis()->SetTitle("Y (mirons)");
histDoseRed->GetZaxis()->SetTitle("Dose (mGy)");
histDoseRed->GetXaxis()->SetLabelSize(0.025);
histDoseRed->GetYaxis()->SetLabelSize(0.025);
histDoseRed->GetZaxis()->SetLabelSize(0.025);
histDoseRed->GetXaxis()->SetTitleSize(0.035);
histDoseRed->GetYaxis()->SetTitleSize(0.035);
histDoseRed->GetZaxis()->SetTitleSize(0.035);
histDoseRed->GetXaxis()->SetTitleOffset(1.4);
histDoseRed->GetYaxis()->SetTitleOffset(1.4);
histDoseRed->GetZaxis()->SetTitleOffset(.6);
histDoseRed->SetTitle("Dose map for RED voxels");

c1->cd(8);
TH2F *histDoseGreen = new TH2F("histDoseGreen","histDoseGreen",100,0,800,100,0,800);
// WARNING : dose scaling to mGy
ntuple2->Draw("y:x>>histDoseGreen","dose/1000","contz");
//gPad->SetLogz();
histDoseGreen->Draw("contz");
histDoseGreen->GetXaxis()->SetTitle("X (microns)");
histDoseGreen->GetYaxis()->SetTitle("Y (mirons)");
histDoseGreen->GetZaxis()->SetTitle("Dose (mGy)");
histDoseGreen->GetXaxis()->SetLabelSize(0.025);
histDoseGreen->GetYaxis()->SetLabelSize(0.025);
histDoseGreen->GetZaxis()->SetLabelSize(0.025);
histDoseGreen->GetXaxis()->SetTitleSize(0.035);
histDoseGreen->GetYaxis()->SetTitleSize(0.035);
histDoseGreen->GetZaxis()->SetTitleSize(0.035);
histDoseGreen->GetXaxis()->SetTitleOffset(1.4);
histDoseGreen->GetYaxis()->SetTitleOffset(1.4);
histDoseGreen->GetZaxis()->SetTitleOffset(.6);
histDoseGreen->SetTitle("Dose map for GREEN voxels");

c1->cd(12);
TH2F *histDoseBlue = new TH2F("histDoseBlue","histDoseBlue",100,0,800,100,0,800);
// WARNING : dose scaling to mGy
ntuple3->Draw("y:x>>histDoseBlue","dose/1000","contz");
//gPad->SetLogz();
histDoseBlue->Draw("contz");
histDoseBlue->GetXaxis()->SetTitle("X (microns)");
histDoseBlue->GetYaxis()->SetTitle("Y (mirons)");
histDoseBlue->GetZaxis()->SetTitle("Dose (mGy)");
histDoseBlue->GetXaxis()->SetLabelSize(0.025);
histDoseBlue->GetYaxis()->SetLabelSize(0.025);
histDoseBlue->GetZaxis()->SetLabelSize(0.025);
histDoseBlue->GetXaxis()->SetTitleSize(0.035);
histDoseBlue->GetYaxis()->SetTitleSize(0.035);
histDoseBlue->GetZaxis()->SetTitleSize(0.035);
histDoseBlue->GetXaxis()->SetTitleOffset(1.4);
histDoseBlue->GetYaxis()->SetTitleOffset(1.4);
histDoseBlue->GetZaxis()->SetTitleOffset(.6);
histDoseBlue->SetTitle("Dose map for BLUE voxels");

}
