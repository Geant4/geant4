// -------------------------------------------------------------------
// $Id: plot.C,v 1.4 2006/06/01 22:25:19 sincerti Exp $
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// This macro needs five files : dose.txt, stoppingPower.txt, range.txt,
// 3DDose.txt and beamPosition.txt
// written by S. Incerti and O. Boissonnade, 10/04/2006
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
	
c1 = new TCanvas ("c1","",20,20,1200,900);
c1.Divide(4,3);

FILE * fp = fopen("dose.txt","r");
Float_t nD,cD;
Int_t ncols=0;
Int_t nlines = 0;

TH1F *h1  = new TH1F("Absorbed dose distribution in Nucleus","Dose distribution in Nucleus",100,0.001,0.5);
TH1F *h10 = new TH1F("Absorbed dose distribution in Cytoplasm","Dose distribution in Cytoplasm",100,0.001,0.2);

while (1) 
   {
      ncols = fscanf(fp,"%f %f",&nD,&cD);
      if (ncols < 0) break;
      h1->Fill(nD);
      h10->Fill(cD);
      nlines++;
   }
fclose(fp);

c1.cd(2);
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
// DOSE IN CYTOPLASM
//*****************

c1.cd(6);
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
    
c1.cd(9);
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
      
c1.cd(10);
  ntuple->Draw("X2:Z2","abs(X2)<50","surf3");
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
  htemp->SetTitle("Range in cell");

//*********************
// INTENSITY HISTOGRAMS 
//*********************

FILE * fp = fopen("phantom.dat","r");
Float_t xVox, yVox, zVox, tmp, den, dose;

Float_t X,Y,Z;
Float_t vox = 0,  mat = 0;
Float_t voxelSizeX, voxelSizeY, voxelSizeZ;

TH1F *h1  = new TH1F("h1","Nucleus marker intensity",100,1,300);
TH1F *h11 = new TH1F("h11 ","",100,1,300);

TH1F *h2  = new TH1F("h2","Cytoplasm marker intensity",100,1,300);
TH1F *h20 = new TH1F("h20 ","",100,1,300);

TNtuple *ntupleYXN = new TNtuple("NUCLEUS","ntuple","Y:X:vox");
TNtuple *ntupleZX = new TNtuple("CYTOPASM","ntuple","Z:X:vox");
TNtuple *ntupleYX = new TNtuple("CYTOPASM","ntuple","Y:X:vox");

nlines=0;
ncols=0;

while (1) 
   {
      if ( nlines == 0 ) ncols = fscanf(fp,"%f %f %f",&tmp,&tmp,&tmp);
      if ( nlines == 1 ) ncols = fscanf(fp,"%f %f %f",&voxelSizeX,&voxelSizeY,&voxelSizeZ);
      if ( nlines == 2 ) ncols = fscanf(fp,"%f %f %f",&tmp,&tmp,&tmp);
      if ( nlines >= 3 ) ncols = fscanf(fp,"%f %f %f %f %f %f", &X, &Y, &Z, &mat, &den, &vox);
      if (ncols < 0) break;

      X= X*voxelSizeX;
      Y= Y*voxelSizeY;
      Z= Z*voxelSizeZ;
  
      if ( mat == 2 )  // noyau
         {
	  if (den==1) h1->Fill( vox );
	  if (den==2) h11->Fill( vox );
	  ntupleYXN->Fill(Y,X,vox);
	 }
      if ( mat == 1 ) // cytoplasm
         {
	  if (den==1) h2->Fill( vox );
	  if (den==2) h20->Fill( vox );
	  ntupleZX->Fill(Z,X,vox);
	  ntupleYX->Fill(Y,X,vox);
	 }
      nlines++;
      
   }
fclose(fp);

// HISTO NUCLEUS

c1.cd(1);
	h1->Draw();
	h1->GetXaxis()->SetLabelSize(0.025);
	h1->GetYaxis()->SetLabelSize(0.025);
	h1->GetXaxis()->SetTitleSize(0.035);
	h1->GetYaxis()->SetTitleSize(0.035);
	h1->GetXaxis()->SetTitleOffset(1.4);
	h1->GetYaxis()->SetTitleOffset(1.4);
	h1->GetXaxis()->SetTitle("Voxel intensity (0-255)");
	h1->GetYaxis()->SetTitle("Number of events");  
	h1->SetLineColor(3);
	h1->SetFillColor(3);   // green

	h11->SetLineColor(8);
	h11->SetFillColor(8);  // dark green
	h11->Draw("same");

// HISTO CYTOPLASM

c1.cd(5);
	h2->Draw();
	h2->GetXaxis()->SetLabelSize(0.025);
	h2->GetYaxis()->SetLabelSize(0.025);
	h2->GetXaxis()->SetTitleSize(0.035);
	h2->GetYaxis()->SetTitleSize(0.035);
	h2->GetXaxis()->SetTitleOffset(1.4);
	h2->GetYaxis()->SetTitleOffset(1.4);
	h2->GetXaxis()->SetTitle("Voxel intensity (0-255)");
	h2->GetYaxis()->SetTitle("Number of events");  
	h2->SetLineColor(2);
	h2->SetFillColor(2);   // red
	
	h20->SetLineColor(5);
	h20->SetFillColor(5);  // yellow (nucleoli)
	h20->Draw("same");

//*************************
// CUMULATED CELL INTENSITY 
//*************************

gStyle->SetOptStat(0000);
gStyle->SetOptFit();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");

//CYTOPLASM

c1.cd(7);  // axe YX
  TH2F *hist = new TH2F("hist","hist",50,-20,20,50,-20,20);
  ntupleYX->Draw("Y:X>>hist","vox","colz");
  gPad->SetLogz();
  hist->Draw("colz");
  hist->GetXaxis()->SetLabelSize(0.025);
  hist->GetYaxis()->SetLabelSize(0.025);
  hist->GetZaxis()->SetLabelSize(0.025);
  hist->GetXaxis()->SetTitleSize(0.035);
  hist->GetYaxis()->SetTitleSize(0.035);
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetTitle("Y (µm)");
  hist->GetYaxis()->SetTitle("X (µm)");
  hist->SetTitle("Cytoplasm intensity on transverse section");

//NUCLEUS

c1.cd(3);  // axe YX
  TH2F *hist = new TH2F("hist","hist",50,-20,20,50,-20,20);
  ntupleYXN->Draw("Y:X>>hist","vox","colz");
  gPad->SetLogz();
  hist->Draw("colz");
  hist->GetXaxis()->SetLabelSize(0.025);
  hist->GetYaxis()->SetLabelSize(0.025);
  hist->GetZaxis()->SetLabelSize(0.025);
  hist->GetXaxis()->SetTitleSize(0.035);
  hist->GetYaxis()->SetTitleSize(0.035);
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetTitle("Y (µm)");
  hist->GetYaxis()->SetTitle("X (µm)");
  hist->SetTitle("Nucleus intensity on transverse section");

//****************
// ENERGY DEPOSITS 
//****************

gStyle->SetOptStat(0000);
gStyle->SetOptFit();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");

FILE * fp = fopen("3DDose.txt","r");

TNtuple *ntuple2 = new TNtuple("CELL","ntuple","yVox:xVox:dose");
TNtuple *ntuple3 = new TNtuple("CELL","ntuple","xVox:zVox:dose");
TNtuple *ntuplezyx = new TNtuple("DOSE","ntuple","zVox:yVox:xVox:dose");

while (1) 
   {
      ncols = fscanf(fp,"%f %f %f %f",&xVox, &yVox, &zVox, &dose);
      if (ncols < 0) break;
      
      xVox= xVox*voxelSizeX;
      yVox= yVox*voxelSizeY;
      zVox= zVox*voxelSizeZ;
      
      ntuple2->Fill(yVox,xVox,dose);
      ntuple3->Fill(xVox,zVox,dose);
      ntuplezyx->Fill(zVox,yVox,xVox,dose);
   }
fclose(fp);

c1.cd(11);
  TH2F *hist = new TH2F("hist","hist",50,-20,20,50,-20,20);
  ntuple2->Draw("yVox:xVox>>hist","dose","contz");
  gPad->SetLogz();
  hist->Draw("contz");
  hist->GetXaxis()->SetLabelSize(0.025);
  hist->GetYaxis()->SetLabelSize(0.025);
  hist->GetZaxis()->SetLabelSize(0.025);
  hist->GetXaxis()->SetTitleSize(0.035);
  hist->GetYaxis()->SetTitleSize(0.035);
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetTitle("Y (µm)");
  hist->GetYaxis()->SetTitle("X (µm)");
  hist->SetTitle("Mean energy deposit -transverse- (z axis in eV)");
  
c1.cd(12);
  TH2F *hist = new TH2F("hist","hist",50,-20,20,50,-20,20);
  ntuple3->Draw("xVox:zVox>>hist","dose","contz");
  gPad->SetLogz();
  hist->Draw("contz");
  hist->GetXaxis()->SetLabelSize(0.025);
  hist->GetYaxis()->SetLabelSize(0.025);
  hist->GetZaxis()->SetLabelSize(0.025);
  hist->GetXaxis()->SetTitleSize(0.035);
  hist->GetYaxis()->SetTitleSize(0.035);
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetTitle("Z (µm)");
  hist->GetYaxis()->SetTitle("X (µm)");
  hist->SetTitle("Mean energy deposit -longitudinal- (z axis in eV)");

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
