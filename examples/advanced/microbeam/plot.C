// -------------------------------------------------------------------
// $Id: plot.C 69157 2013-04-19 14:05:11Z gcosmo $
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
c1->Divide(4,3);

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
TNtuple *ntupleZX = new TNtuple("CYTOPLASM","ntuple","Z:X:vox");
TNtuple *ntupleYX = new TNtuple("CYTOPLASM","ntuple","Y:X:vox");

Int_t nlines=0;
Int_t ncols=0;

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

c1->cd(1);
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

c1->cd(5);
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

c1->cd(7);  // axe YX
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
  hist->GetXaxis()->SetTitle("Y (um)");
  hist->GetYaxis()->SetTitle("X (um)");
  hist->SetTitle("Cytoplasm intensity on transverse section");

//NUCLEUS

c1->cd(3);  // axe YX
  TH2F *hist2 = new TH2F("hist2","hist2",50,-20,20,50,-20,20);
  ntupleYXN->Draw("Y:X>>hist2","vox","colz");
  gPad->SetLogz();
  hist2->Draw("colz");
  hist2->GetXaxis()->SetLabelSize(0.025);
  hist2->GetYaxis()->SetLabelSize(0.025);
  hist2->GetZaxis()->SetLabelSize(0.025);
  hist2->GetXaxis()->SetTitleSize(0.035);
  hist2->GetYaxis()->SetTitleSize(0.035);
  hist2->GetXaxis()->SetTitleOffset(1.4);
  hist2->GetYaxis()->SetTitleOffset(1.4);
  hist2->GetXaxis()->SetTitle("Y (um)");
  hist2->GetYaxis()->SetTitle("X (um)");
  hist2->SetTitle("Nucleus intensity on transverse section");

//

system ("rm -rf microbeam.root");
system ("hadd -O microbeam.root microbeam_*.root");

TFile f("microbeam.root"); 

TNtuple* ntuple0;
TNtuple* ntuple1;
TNtuple* ntuple2;
TNtuple* ntuple3;
TNtuple* ntuple4;

ntuple0 = (TNtuple*)f.Get("ntuple0"); 
ntuple1 = (TNtuple*)f.Get("ntuple1"); 
ntuple2 = (TNtuple*)f.Get("ntuple2"); 
ntuple3 = (TNtuple*)f.Get("ntuple3"); 
ntuple4 = (TNtuple*)f.Get("ntuple4"); 

TH1F *h1bis  = new TH1F("h1bis","Dose distribution in Nucleus",100,0.001,1.);
TH1F *h10 = new TH1F("h10bis","Dose distribution in Cytoplasm",100,0.001,.2);

c1->cd(2);

        ntuple3->Project("h1bis","doseN");
	scale = 1/h1bis->Integral();
	h1bis->Scale(scale);
	h1bis->Draw();
	h1bis->GetXaxis()->SetLabelSize(0.025);
	h1bis->GetYaxis()->SetLabelSize(0.025);
	h1bis->GetXaxis()->SetTitleSize(0.035);
	h1bis->GetYaxis()->SetTitleSize(0.035);
	h1bis->GetXaxis()->SetTitleOffset(1.4);
	h1bis->GetYaxis()->SetTitleOffset(1.4);
	h1bis->GetXaxis()->SetTitle("Absorbed dose (Gy)");
	h1bis->GetYaxis()->SetTitle("Fraction of events");
	h1bis->SetLineColor(3);
	h1bis->SetFillColor(3);

//*****************
// DOSE IN CYTOPLASM
//*****************

c1->cd(6);
        ntuple3->Project("h10bis","doseC");
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

TH1F *h2bis = new TH1F("h2bis","Beam stopping power at cell entrance",200,0,300); 
    
c1->cd(9);
        ntuple0->Project("h2bis","sp");
	scale = 1/h2bis->Integral();
	h2bis->Scale(scale);
	h2bis->Draw();
	h2bis->GetXaxis()->SetLabelSize(0.025);
	h2bis->GetYaxis()->SetLabelSize(0.025);
	h2bis->GetXaxis()->SetTitleSize(0.035);
	h2bis->GetYaxis()->SetTitleSize(0.035);
	h2bis->GetXaxis()->SetTitleOffset(1.4);
	h2bis->GetYaxis()->SetTitleOffset(1.4);
  	h2bis->GetXaxis()->SetTitle("dE/dx (keV/um)");
	h2bis->GetYaxis()->SetTitle("Fraction of events");
  	h2bis->SetTitle("dE/dx at cell entrance");
	h2bis->SetFillColor(4);
	h2bis->SetLineColor(4);
	h2bis->Fit("gaus");
	gaus->SetLineColor(6);
	h2bis->Fit("gaus");

//**************
// RANGE IN CELL
//**************

Double_t Xc,Zc,X1,Y1,Z1,X2,Y2,Z2;

// X position of target in World
Xc = -1295.59e3 - 955e3*sin(10*TMath::Pi()/180); 

// Z position of target in World
Zc = -1327e3 + 955e3*cos(10*TMath::Pi()/180); 

// Line alignment (cf MicrobeamEMField.cc)
Xc = Xc + 5.24*cos(10*TMath::Pi()/180);
Zc = Zc + 5.24*sin(10*TMath::Pi()/180);

TNtuple *ntupleR = new TNtuple("Rmax","ntuple","Z2:Y2:X2");
Double_t x,y,z,xx,zz;
ntuple2->SetBranchAddress("x",&x);
ntuple2->SetBranchAddress("y",&y);
ntuple2->SetBranchAddress("z",&z);
Int_t nentries = (Int_t)ntuple2->GetEntries();
for (Int_t i=0;i<nentries;i++) 
{
      ntuple2->GetEntry(i);
      X1=x;
      Y1=y;
      Z1=z;
      xx = X1-Xc;
      zz = Z1-Zc;
      Z2 = zz*cos(10*TMath::Pi()/180)-xx*sin(10*TMath::Pi()/180);
      X2 = zz*sin(10*TMath::Pi()/180)+xx*cos(10*TMath::Pi()/180);
      Y2 = Y1;    
      ntupleR->Fill(Z2,Y2,X2);
}
     
c1->cd(10);
  ntupleR->Draw("X2:Z2","abs(X2)<50","surf3");
  gPad->SetLogz();

//****************
// ENERGY DEPOSITS 
//****************

gStyle->SetOptStat(0000);
gStyle->SetOptFit();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");

c1->cd(11);
  TH2F *histbis = new TH2F("histbis","histbis",50,-20,20,50,-20,20);
  ntuple4->Draw("y*0.359060:x*0.359060>>histbis","doseV","contz");
  gPad->SetLogz();
  histbis->Draw("contz");
  histbis->GetXaxis()->SetLabelSize(0.025);
  histbis->GetYaxis()->SetLabelSize(0.025);
  histbis->GetZaxis()->SetLabelSize(0.025);
  histbis->GetXaxis()->SetTitleSize(0.035);
  histbis->GetYaxis()->SetTitleSize(0.035);
  histbis->GetXaxis()->SetTitleOffset(1.4);
  histbis->GetYaxis()->SetTitleOffset(1.4);
  histbis->GetXaxis()->SetTitle("Y (um)");
  histbis->GetYaxis()->SetTitle("X (um)");
  histbis->SetTitle("Mean energy deposit -transverse- (z axis in eV)");

c1->cd(12);
  TH2F *histter = new TH2F("histter","histter",50,-20,20,50,-20,20);
  ntuple4->Draw("x*0.359060:(z+1500/0.162810+21)*0.162810>>histter","doseV","contz");
  gPad->SetLogz();
  histter->Draw("contz");
  histter->GetXaxis()->SetLabelSize(0.025);
  histter->GetYaxis()->SetLabelSize(0.025);
  histter->GetZaxis()->SetLabelSize(0.025);
  histter->GetXaxis()->SetTitleSize(0.035);
  histter->GetYaxis()->SetTitleSize(0.035);
  histter->GetXaxis()->SetTitleOffset(1.4);
  histter->GetYaxis()->SetTitleOffset(1.4);
  histter->GetXaxis()->SetTitle("Z (um)");
  histter->GetYaxis()->SetTitle("X (um)");
  histter->SetTitle("Mean energy deposit -longitudinal- (z axis in eV)");

//*******************************
// BEAM POSITION AT CELL ENTRANCE
//*******************************
 
gStyle->SetOptStat(0000);
gStyle->SetOptFit();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");

TH1F *h77 = new TH1F("hx","h1",200,-10,10); 
TH1F *h88 = new TH1F("hy","h1",200,-10,10); 
   
c1->cd(4);
        ntuple1->Project("hx","x");
	scale = 1/h77->Integral();
	h77->Scale(scale);
	h77->Draw();
	h77->GetXaxis()->SetLabelSize(0.025);
	h77->GetYaxis()->SetLabelSize(0.025);
	h77->GetXaxis()->SetTitleSize(0.035);
	h77->GetYaxis()->SetTitleSize(0.035);
	h77->GetXaxis()->SetTitleOffset(1.4);
	h77->GetYaxis()->SetTitleOffset(1.4);
  	h77->GetXaxis()->SetTitle("Position (um)");
	h77->GetYaxis()->SetTitle("Fraction of events");
  	h77->SetTitle("Beam X position on cell");
	h77->SetFillColor(4);
	h77->SetLineColor(4);
	h77->Fit("gaus");

c1->cd(8);
        ntuple1->Project("hy","y");
	scale = 1/h88->Integral();
	h88->Scale(scale);
	h88->Draw();
	h88->GetXaxis()->SetLabelSize(0.025);
	h88->GetYaxis()->SetLabelSize(0.025);
	h88->GetXaxis()->SetTitleSize(0.035);
	h88->GetYaxis()->SetTitleSize(0.035);
	h88->GetXaxis()->SetTitleOffset(1.4);
	h88->GetYaxis()->SetTitleOffset(1.4);
  	h88->GetXaxis()->SetTitle("Position (um)");
	h88->GetYaxis()->SetTitle("Fraction of events");
  	h88->SetTitle("Beam Y position on cell");
	h88->SetFillColor(4);
	h88->SetLineColor(4);
	h88->Fit("gaus");
}
