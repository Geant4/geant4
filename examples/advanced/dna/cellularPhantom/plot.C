// -------------------------------------------------------------------
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT,
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// Written by S. Incerti, 10/09/2024
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

// IF no merging active in simulation
//system ("rm -rf phantom.root");
//system ("hadd -O phantom.root phantom_t*.root");

TFile *f = new TFile ("phantom.root");

TNtuple* ntuple1;
TNtuple* ntuple2;
TNtuple* ntuple3;

ntuple1 = (TNtuple*)f->Get("ntuple1");
ntuple2 = (TNtuple*)f->Get("ntuple2");
ntuple3 = (TNtuple*)f->Get("ntuple3");

//----------------------
// Sum of ntuples
//----------------------

Double_t * tabVoxelXRed = new Double_t [numberVoxTot];
Double_t * tabVoxelXGreen = new Double_t [numberVoxTot];
Double_t * tabVoxelXBlue = new Double_t [numberVoxTot];

Double_t * tabVoxelYRed = new Double_t [numberVoxTot];
Double_t * tabVoxelYGreen = new Double_t [numberVoxTot];
Double_t * tabVoxelYBlue = new Double_t [numberVoxTot];

Double_t * tabVoxelZRed = new Double_t [numberVoxTot];
Double_t * tabVoxelZGreen = new Double_t [numberVoxTot];
Double_t * tabVoxelZBlue = new Double_t [numberVoxTot];

Double_t * tabVoxelEnergyRed = new Double_t [numberVoxTot];
Double_t * tabVoxelEnergyGreen = new Double_t [numberVoxTot];
Double_t * tabVoxelEnergyBlue = new Double_t [numberVoxTot];

Double_t * tabVoxelDoseRed = new Double_t [numberVoxTot];
Double_t * tabVoxelDoseGreen = new Double_t [numberVoxTot];
Double_t * tabVoxelDoseBlue = new Double_t [numberVoxTot];

// Initialisation of the arrays
for (Int_t i = 0; i < numberVoxRed; i++)
{
  tabVoxelXRed[i] = 0;
  tabVoxelYRed[i] = 0;
  tabVoxelZRed[i] = 0;
  tabVoxelEnergyRed[i] = 0;
  tabVoxelDoseRed[i] = 0;
}
for (Int_t i = 0; i < numberVoxGreen; i++)
{
  tabVoxelXGreen[i] = 0;
  tabVoxelYGreen[i] = 0;
  tabVoxelZGreen[i] = 0;
  tabVoxelEnergyGreen[i] = 0;
  tabVoxelDoseGreen[i] = 0;
}
for (Int_t i = 0; i < numberVoxBlue; i++)
{
  tabVoxelXBlue[i] = 0;
  tabVoxelYBlue[i] = 0;
  tabVoxelZBlue[i] = 0;
  tabVoxelEnergyBlue[i] = 0;
  tabVoxelDoseBlue[i] = 0;
}

Double_t x, y, z, energy, dose;
Int_t voxelID;
Double_t nrjRed=0.;
Double_t nrjGreen=0.;
Double_t nrjBlue=0.;
Double_t doseRed=0.;
Double_t doseGreen=0.;
Double_t doseBlue=0.;

//

ntuple1->SetBranchAddress("x",&x);
ntuple1->SetBranchAddress("y",&y);
ntuple1->SetBranchAddress("z",&z);
ntuple1->SetBranchAddress("energy",&energy);
ntuple1->SetBranchAddress("dose",&dose);
ntuple1->SetBranchAddress("voxelID",&voxelID);

// RED

Long_t nentriesRed = (Long_t)ntuple1->GetEntries();
for (Long_t i=0;i<nentriesRed;i++)
{
      x=0;
      y=0;
      z=0;
      energy=0;
      dose=0;
      voxelID=0;

      ntuple1->GetEntry(i);
      if (energy > 0)
      {
        nrjRed=nrjRed+energy;
        doseRed=doseRed+dose;

        tabVoxelXRed[voxelID] = x;
        tabVoxelYRed[voxelID] = y;
        tabVoxelZRed[voxelID] = z;
        tabVoxelEnergyRed[voxelID] = tabVoxelEnergyRed[voxelID] + energy;
        tabVoxelDoseRed[voxelID] = tabVoxelDoseRed[voxelID] + dose;
      }
}

ntuple2->SetBranchAddress("x",&x);
ntuple2->SetBranchAddress("y",&y);
ntuple2->SetBranchAddress("z",&z);
ntuple2->SetBranchAddress("energy",&energy);
ntuple2->SetBranchAddress("dose",&dose);
ntuple2->SetBranchAddress("voxelID",&voxelID);

// GREEN

Long_t nentriesGreen = (Long_t)ntuple2->GetEntries();
for (Long_t i=0;i<nentriesGreen;i++)
{
      x=0;
      y=0;
      z=0;
      energy=0;
      dose=0;
      voxelID=0;

      ntuple2->GetEntry(i);
      if (energy > 0)
      {
        nrjGreen=nrjGreen+energy;
        doseGreen=doseGreen+dose;

        tabVoxelXGreen[voxelID] = x;
        tabVoxelYGreen[voxelID] = y;
        tabVoxelZGreen[voxelID] = z;
        tabVoxelEnergyGreen[voxelID] = tabVoxelEnergyGreen[voxelID] + energy;
        tabVoxelDoseGreen[voxelID] = tabVoxelDoseGreen[voxelID] + dose;
      }
}

// BLUE

ntuple3->SetBranchAddress("x",&x);
ntuple3->SetBranchAddress("y",&y);
ntuple3->SetBranchAddress("z",&z);
ntuple3->SetBranchAddress("energy",&energy);
ntuple3->SetBranchAddress("dose",&dose);
ntuple3->SetBranchAddress("voxelID",&voxelID);

Long_t nentriesBlue = (Long_t)ntuple3->GetEntries();
for (Long_t i=0;i<nentriesBlue;i++)
{
      x=0;
      y=0;
      z=0;
      energy=0;
      dose=0;
      voxelID=0;

      ntuple3->GetEntry(i);
      if (energy > 0)
      {
        nrjBlue=nrjBlue+energy;
        doseBlue=doseBlue+dose;
        tabVoxelXBlue[voxelID] = x;
        tabVoxelYBlue[voxelID] = y;
        tabVoxelZBlue[voxelID] = z;
        tabVoxelEnergyBlue[voxelID] = tabVoxelEnergyBlue[voxelID] + energy;
        tabVoxelDoseBlue[voxelID] = tabVoxelDoseBlue[voxelID] + dose;
      }
}

// To liberate memory
f->Close();

TFile *f2 = new TFile ("results.root","RECREATE");
//

TNtuple *ntupleRED = new TNtuple ("RED","RED","x:y:z:energy:dose");
TNtuple *ntupleGREEN = new TNtuple ("GREEN","GREEN","x:y:z:energy:dose");
TNtuple *ntupleBLUE = new TNtuple ("BLUE","BLUE","x:y:z:energy:dose");

// Global sums
for (Int_t i = 0; i < numberVoxTot; i++)
{
  ntupleRED->Fill(tabVoxelXRed[i],tabVoxelYRed[i],tabVoxelZRed[i],tabVoxelEnergyRed[i],tabVoxelDoseRed[i]);
}
for (Int_t i = 0; i < numberVoxTot; i++)
{
  ntupleGREEN->Fill(tabVoxelXGreen[i],tabVoxelYGreen[i],tabVoxelZGreen[i],tabVoxelEnergyGreen[i],tabVoxelDoseGreen[i]);
}
for (Int_t i = 0; i < numberVoxTot; i++)
{
  ntupleBLUE->Fill(tabVoxelXBlue[i],tabVoxelYBlue[i],tabVoxelZBlue[i],tabVoxelEnergyBlue[i],tabVoxelDoseBlue[i]);
}

//---------------------------------
// Absorbed energy distributions
//---------------------------------

c1->cd(2);
gPad->SetLogy();
ntupleRED->Draw("energy","energy>0");
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
ntupleGREEN->Draw("energy","energy>0");
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
ntupleBLUE->Draw("energy","energy>0");
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
ntupleRED->Draw("y:x>>histNrjRed","energy","contz");
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
ntupleGREEN->Draw("y:x>>histNrjGreen","energy","contz");
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
ntupleBLUE->Draw("y:x>>histNrjBlue","energy","contz");
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
ntupleRED->Draw("y:x>>histDoseRed","dose/1000","contz");
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
ntupleGREEN->Draw("y:x>>histDoseGreen","dose/1000","contz");
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
ntupleBLUE->Draw("y:x>>histDoseBlue","dose/1000","contz");
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

//----------------------------
// SUMMARY
//----------------------------

cout << endl;
cout << "- Summary --------------------------------------------------" << endl;
cout << endl;
cout << "  Total number of voxels in phantom           = " << numberVoxTot << endl;
cout << "  Total number of RED voxels in phantom       = " << numberVoxRed << endl;
cout << "  Total number of GREEN voxels in phantom     = " << numberVoxGreen << endl;
cout << "  Total number of BLUE voxels in phantom      = " << numberVoxBlue << endl;
cout << endl;
cout << "  Total absorbed energy in RED   voxels (MeV) = "  << nrjRed/1E3 << endl;
cout << "  Total absorbed energy in GREEN voxels (MeV) = "  << nrjGreen/1E3 << endl;
cout << "  Total absorbed energy in BLUE  voxels (MeV) = "  << nrjBlue/1E3 << endl;
cout << endl;
cout << "  Total absorbed dose   in RED   voxels (Gy)  = "  << doseRed << endl;
cout << "  Total absorbed dose   in GREEN voxels (Gy)  = "  << doseGreen << endl;
cout << "  Total absorbed dose   in BLUE  voxels (Gy)  = "  << doseBlue << endl;
cout << endl;
cout << "------------------------------------------------------------" << endl;

// End
f2->Write();

}
