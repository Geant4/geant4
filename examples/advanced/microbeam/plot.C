// -------------------------------------------------------------------
// $Id: plot.C,v 1.2 2006-04-06 15:22:36 sincerti Exp $
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// This macro needs five files : dose.txt, stoppingPower.txt, range.txt,
// 3DDose.txt and beamPosition.txt
// written by S. Incerti, 01/04/2006
// *********************************************************************
{
gROOT->Reset();

//****************
// DOSE IN NUCLEUS
//****************

gStyle->SetOptStat(0000);
gStyle->SetOptFit();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
Double_t scale;
	
c1 = new TCanvas ("c1","",200,20,1200,600);
c1.Divide(4,2);

FILE * fp = fopen("dose.txt","r");
Float_t nD,cD;
Int_t ncols;
Int_t nlines = 0;
TH1F *h1  = new TH1F("Dose distribution in Nucleus","Dose distribution in Nucleus",100,0.001,0.5);
TH1F *h10 = new TH1F("Dose distribution in Cytoplasm","Dose distribution in Cytoplasm",100,0.001,0.5);

while (1) 
   {
      ncols = fscanf(fp,"%f %f",&nD,&cD);
      if (ncols < 0) break;
      h1->Fill(nD);
      h10->Fill(cD);
      nlines++;
   }
fclose(fp);

c1.cd(1);
        scale = 1/h1->Integral();
	h1->Scale(scale);
	h1->Draw();
	h1->GetXaxis()->SetLabelSize(0.025);
	h1->GetYaxis()->SetLabelSize(0.025);
	h1->GetXaxis()->SetTitleSize(0.035);
	h1->GetYaxis()->SetTitleSize(0.035);
	h1->GetXaxis()->SetTitleOffset(1.4);
	h1->GetYaxis()->SetTitleOffset(1.4);
	h1->GetXaxis()->SetTitle("Absorbed dose (Gy)");
	h1->GetYaxis()->SetTitle("Fraction of events");
	h1->SetLineColor(3);
	h1->SetFillColor(3);

//*****************
// DOSE IN CTOPLASM
//*****************

c1.cd(5);
        scale = 1/h10->Integral();
	h10->Scale(scale);
	h10->Draw();
	h10->GetXaxis()->SetLabelSize(0.025);
	h10->GetYaxis()->SetLabelSize(0.025);
	h10->GetXaxis()->SetTitleSize(0.035);
	h10->GetYaxis()->SetTitleSize(0.035);
	h10->GetXaxis()->SetTitleOffset(1.4);
	h10->GetYaxis()->SetTitleOffset(1.4);
	h10->GetXaxis()->SetTitle("Absorbed dose (Gy)");
	h10->GetYaxis()->SetTitle("Fraction of events");
	h10->SetLineColor(2);
	h10->SetFillColor(2);

//********************************
// STOPPING POWER AT CELL ENTRANCE
//********************************
 
gStyle->SetOptStat(0000);
gStyle->SetOptFit();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");

Float_t d;
FILE * fp = fopen("stoppingPower.txt","r");

TH1F *h2 = new TH1F("Beam stopping Power at cell entrance","h1",200,0,300); 
while (1) 
   {
      ncols = fscanf(fp,"%f",&d);
      if (ncols < 0) break;
      h2->Fill(d);
      nlines++;
    }
fclose(fp);
    
c1.cd(2);
        scale = 1/h2->Integral();
	h2->Scale(scale);
	h2->Draw();
	h2->GetXaxis()->SetLabelSize(0.025);
	h2->GetYaxis()->SetLabelSize(0.025);
	h2->GetXaxis()->SetTitleSize(0.035);
	h2->GetYaxis()->SetTitleSize(0.035);
	h2->GetXaxis()->SetTitleOffset(1.4);
	h2->GetYaxis()->SetTitleOffset(1.4);
  	h2->GetXaxis()->SetTitle("dE/dx (keV/µm)");
	h2->GetYaxis()->SetTitle("Fraction of events");
  	h2->SetTitle("dE/dx at cell entrance");
	h2->SetFillColor(4);
	h2->SetLineColor(4);
	gaus->SetLineColor(6);
	h2->Fit("gaus");

//**************
// RANGE IN CELL
//**************

Float_t Xc,Zc,X1,Y1,Z1,x,z,X2,Y2,Z2;
Float_t d;

// X position of target in World
Xc = -1295.59e3 - 955e3*sin(10*TMath::Pi()/180); 

// Z position of target in World
Zc = -1327e3 + 955e3*cos(10*TMath::Pi()/180); 

// Line alignment (cf MicrobeamEMField.cc)
Xc = Xc + 5.24*cos(10*TMath::Pi()/180);
Zc = Zc + 5.24*sin(10*TMath::Pi()/180);

FILE * fp = fopen("range.txt","r");

TNtuple *ntuple = new TNtuple("Rmax","ntuple","Z2:Y2:X2");

while (1) 
   {
      ncols = fscanf(fp,"%f %f %f",&X1,&Y1,&Z1);
      if (ncols < 0) break;
      x = X1-Xc;
      z = Z1-Zc;
      Z2 = z*cos(10*TMath::Pi()/180)-x*sin(10*TMath::Pi()/180);
      X2 = z*sin(10*TMath::Pi()/180)+x*cos(10*TMath::Pi()/180);
      Y2 = Y1;
      
      ntuple->Fill(Z2,Y2,X2);
      nlines++;
   }
fclose(fp);
      
c1.cd(6);
  ntuple->Draw("X2:Z2","abs(X2)<50","colz");
  gPad->SetLogz();
  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetZaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.4);
  htemp->GetYaxis()->SetTitleOffset(1.4);
  htemp->GetXaxis()->SetTitle("Z (µm)");
  htemp->GetYaxis()->SetTitle("X (µm)");
  htemp->SetTitle("Range in cell phantom");

//*****************************
// ENERGY DEPOSITS - HORIZONTAL
//*****************************

Float_t xVox, yVox, zVox, dose, tmp;

Float_t voxelSizeX, voxelSizeY, voxelSizeZ;
FILE * fp = fopen("phantom.dat","r");
ncols = fscanf(fp,"%f %f %f",&tmp,&tmp,&tmp);
ncols = fscanf(fp,"%f %f %f",&voxelSizeX,&voxelSizeY,&voxelSizeZ);
fclose(fp);

FILE * fp = fopen("3DDose.txt","r");

TNtuple *ntupleDose = new TNtuple("Dose","ntuple","xVox:yVox:zVox:dose");

TH2F *hxy = new TH2F("xy","xy",128,-64,63,128,-64,63);  
TH2F *hxz = new TH2F("xz","xz",128,-64,63,59,-29,28);  
TH2F *hyz = new TH2F("yz","yz",128,-64,63,59,-29,28);  

while (1) 
   {
      ncols = fscanf(fp,"%f %f %f %f",&xVox,&yVox,&zVox,&dose);
      if (ncols < 0) break;

      voxelSizeX=1;
      voxelSizeY=1;
      voxelSizeZ=1;

      xVox = xVox * voxelSizeX;
      yVox = yVox * voxelSizeY;
      zVox = zVox * voxelSizeZ;

      ntupleDose->Fill(xVox,yVox,zVox,dose);
      hxy->Fill(xVox,yVox);
      
      nlines++;
   }
fclose(fp);
     
c1.cd(7);
  hxy->Draw("colz");
  gPad->SetLogz();

  hxy->GetXaxis()->SetLabelSize(0.025);
  hxy->GetYaxis()->SetLabelSize(0.025);
  hxy->GetZaxis()->SetLabelSize(0.025);
  hxy->GetXaxis()->SetTitleSize(0.035);
  hxy->GetYaxis()->SetTitleSize(0.035);
  hxy->GetXaxis()->SetTitleOffset(1.4);
  hxy->GetYaxis()->SetTitleOffset(1.4);
  hxy->GetXaxis()->SetTitle("Y (voxels)");
  hxy->GetYaxis()->SetTitle("X (voxels)");
  hxy->SetTitle("Cell density cumulated on transverse section");


c1.cd(3);

  gStyle->SetPalette(1);
  ntupleDose->Draw("dose:xVox:yVox","dose>0","iso");
  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetZaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetZaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.4);
  htemp->GetYaxis()->SetTitleOffset(1.4);
  htemp->GetZaxis()->SetTitleOffset(1.4);
  htemp->GetXaxis()->SetTitle("Y (voxels)");
  htemp->GetYaxis()->SetTitle("X (voxels)");
  htemp->GetZaxis()->SetTitle("Energy deposit (eV)");
  htemp->SetTitle("Average energy deposit per voxel cumulated on cell transverse section");

/*
  ntupleDose->Draw("dose>>htemp","dose>0","");
  gPad->SetLogz();
  scale = 1/htemp->Integral();
  htemp->Scale(scale);
  htemp->Draw();
  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetZaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.4);
  htemp->GetYaxis()->SetTitleOffset(1.4);
  htemp->GetXaxis()->SetTitle("Energy deposit (eV)");
  htemp->GetYaxis()->SetTitle("Fraction of events");
  htemp->SetTitle("Average energy deposit in voxels");
  htemp->SetFillColor(6);
  htemp->SetLineColor(6);
  gPad->SetLogy();
  */
  
//*******************************
// BEAM POSITION AT CELL ENTRANCE
//*******************************
 
gStyle->SetOptStat(0000);
gStyle->SetOptFit();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");

Float_t bx, by;
FILE * fp = fopen("beamPosition.txt","r");

TH1F *h77 = new TH1F("Beam X transverse position at cell entrance","h1",200,-10,10); 
TH1F *h88 = new TH1F("Beam Y transverse position at cell entrance","h1",200,-10,10); 
while (1) 
   {
      ncols = fscanf(fp,"%f %f",&bx, &by);
      if (ncols < 0) break;
      h77->Fill(bx);
      h88->Fill(by);
      nlines++;
    }
fclose(fp);
    
c1.cd(4);
        scale = 1/h77->Integral();
	h77->Scale(scale);
	h77->Draw();
	h77->GetXaxis()->SetLabelSize(0.025);
	h77->GetYaxis()->SetLabelSize(0.025);
	h77->GetXaxis()->SetTitleSize(0.035);
	h77->GetYaxis()->SetTitleSize(0.035);
	h77->GetXaxis()->SetTitleOffset(1.4);
	h77->GetYaxis()->SetTitleOffset(1.4);
  	h77->GetXaxis()->SetTitle("Position (µm)");
	h77->GetYaxis()->SetTitle("Fraction of events");
  	h77->SetTitle("Beam X position on cell");
	h77->SetFillColor(4);
	h77->SetLineColor(4);
	gaus->SetLineColor(6);
	h77->Fit("gaus");

c1.cd(8);
        scale = 1/h88->Integral();
	h88->Scale(scale);
	h88->Draw();
	h88->GetXaxis()->SetLabelSize(0.025);
	h88->GetYaxis()->SetLabelSize(0.025);
	h88->GetXaxis()->SetTitleSize(0.035);
	h88->GetYaxis()->SetTitleSize(0.035);
	h88->GetXaxis()->SetTitleOffset(1.4);
	h88->GetYaxis()->SetTitleOffset(1.4);
  	h88->GetXaxis()->SetTitle("Position (µm)");
	h88->GetYaxis()->SetTitle("Fraction of events");
  	h88->SetTitle("Beam Y position on cell");
	h88->SetFillColor(4);
	h88->SetLineColor(4);
	gaus->SetLineColor(6);
	h88->Fit("gaus");



}
