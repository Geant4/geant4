// Requires: test/tmp/run2.root datafile (generate by ./run.sh run2)
// Original version: Alain Boudard
// Modified by: Pekka Kaitaniemi 23.1.2008: code cleanup for Geant4 release

// SATURNE experiment (n double diff cross at various angles)
//  A.B. 8/10/2007

//#ifndef _ISA_TROLL_	/* includes specific to isa_troll */
//#define _ISA_TROLL_

//#include <stdio.h>
//#define 99 99
//#include <string.h>
//#include <stdlib.h>
//#include <ctype.h>
//#include <math.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "Riostream.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TCut.h"
#include "TTimeStamp.h"
#include "TPaveText.h"
#include "TPaveLabel.h"
#include "TText.h"
#include "TCanvas.h"
//#endif

void protonddrun1() 
{
  gROOT->ProcessLine(".x scripts/rootlogon.C");
  gROOT->SetStyle("clearRetro");
  nCrossSection();
}

void plotExpTheta(Char_t* file_name, Char_t* racine, Float_t fnor)
  // Read experimental data from file: file_name (in path: racine)

  // Cross sections and errors are multiplied by fnor and plotted in a
  // TGraphErrors
{
  // Experimental points
  Char_t toto[256];
  Int_t e;
  Float_t e_n[100];
  Float_t sig[100];
  Float_t dsig[100];
  Float_t de_n[100];
  Int_t c,iexp=0;

  //Concatenation Racine+File_Name
  e=sprintf(toto,"%s%s",racine,file_name);
  cout<<toto<<endl;

  FILE* fexp160 =  fopen(toto,"r");

  iexp=0;
  do {
    c=fscanf(fexp160,"%f %f %f %f ",&e_n[iexp],&de_n[iexp],&sig[iexp],&dsig[iexp]);
    if(c==EOF) break;
    sig[iexp]=fnor*sig[iexp];
    dsig[iexp]=fnor*dsig[iexp];
    //	cout<<e_n[iexp]<<"  "<<sig[iexp]<<"  "<<dsig[iexp]<<endl;
    iexp=iexp+1;
    if(iexp>100) cout<<"More experimental values than the dimension"<<endl;
  }
  while(1);

  TGraphErrors* grexp160= new TGraphErrors(iexp,e_n,sig,de_n,dsig);
  grexp160->SetMarkerColor(4);
  grexp160->SetMarkerStyle(21);
  grexp160->SetMarkerSize(0.4);

  grexp160->Draw("PZ");
  return;
}

void plotTheorTheta(TTree* ref,Char_t* titre,TCanvas* c1,Double_t emin,Double_t emax,Int_t logE,TCut neutron,TCut thet0,Double_t fnor0,Double_t fnora,Char_t* angle,Int_t first, Int_t color = kRed)
  // Project histos (TH1F) using variables from a TTree structure.
  // 100 channels from emin to emax logE =1 for log scale in X (0 for
  // linear scale) neutron and thet0 are 2 cuts applied A
  // normalisation fnor0*fnora is applied (fnora is used for the
  // legend (angle) positionning) first=1 is used for the first histo
  // (proper use of "same" for others)
{
  fnor0=fnor0*fnora;
  Double_t de;
  TH1F* hist0;
  TH1F* hist_de;

  if (logE) { // Logarithmic plot in X
    c1->SetLogx(); 
    Float_t xbins[99+1];
    Float_t fact=(log(emax)-log(emin))/99;
    for(Int_t i=0;i<99+1;++i) {xbins[i]=exp(log(emin)+fact*i);}	
    hist0 = new TH1F("hist0","",99,xbins);
    hist_de = new TH1F("hist_de","",99,xbins);
    for(Int_t i=0;i<99;++i) {hist_de->SetBinContent(i,(Double_t)1./(xbins[i+1]-xbins[i]));}
    de=1.0;
  } else { // Linear plot in X
    hist0 = new TH1F("hist0","",99,emin,emax);
    de=(emax-emin)/99;
  }

  hist0->SetFillStyle(0);
  hist0->SetStats(kFALSE);
  hist0->SetTitle(titre);
  hist0->GetXaxis()->SetTitle("Proton Energy (MeV)");
  hist0->GetXaxis()->CenterTitle(true);
  hist0->GetXaxis()->SetLabelSize(0.03);
  hist0->GetYaxis()->SetTitle("Cross section (mb)");
  hist0->GetYaxis()->CenterTitle(true);
  hist0->GetYaxis()->SetLabelSize(0.03);
  hist0->GetYaxis()->SetTitleOffset(1.3);
  hist0->SetLineColor(color);
	
  // Project in the hist histo
  if(first){ref->Draw("Enerj >> hist0",neutron && thet0);}
  else {ref->Draw("Enerj >> hist0",neutron && thet0,"same");}
	
  fnor0=fnor0/de;
  hist0->Scale(fnor0);
  if (logE) {hist0->Multiply(hist0,hist_de,1.,1.);	
  // Full y extension of the picture:
  hist0->SetMinimum(1.e-13);
  hist0->SetMaximum(1.e+3);
  // Legend (angle):
  Float_t y1=200.0 + 40.0;
  Float_t y2=400.0 + 40.0;
  y1=y1*fnora;
  y2=y2*fnora;	
  TPaveLabel* p1 = new TPaveLabel(1.5,y1,2.5,y2,angle);
  p1->SetBorderSize(0);
  p1->SetFillColor(0);
  p1->SetTextSize(2.0);
  p1->Draw("same");
  delete hist_de;}
  else {
    // Full y extension of the picture:
    hist0->SetMinimum(1.e-3);
    hist0->SetMaximum(1.e+1);
    // Legend (angle):
    Float_t y1=3.;
    Float_t y2=4.;
    y1=y1*fnora;
    y2=y2*fnora;	
    TPaveLabel* p1 = new TPaveLabel(200,y1,300,y2,angle);
    p1->SetBorderSize(0);
    p1->SetFillColor(0);
    p1->SetTextSize(0.8);
    p1->Draw("same");}
	
  //	hist0->DrawCopy();
  return;
}

void nCrossSection(){
  TFile* fref = new TFile("tmp/run1.root");
  TFile* ffref = new TFile("tmp/run1ref.root");

  Char_t* titre="p(1.2 GeV) + 208Pb (INCL4+ABLA)";
	
  Char_t* psFileName="p_cross_section_run1.ps";

  Char_t* racine="./data/proton/pb/"; // Path to experimental files

  // Read the ROOT tree to different variables so that we can handle
  // them easily.
  TTree* ref = (TTree*) fref->Get("h101");
  TTree* reffort = (TTree*) ffref->Get("h101");
	
  // Prepare the postscript file name for several pages:
  Char_t psFileName_debut[256];
  Char_t psFileName_fin[256];
  Int_t e;
  Char_t* debut ="(";
  Char_t* fin =")";
  e=sprintf(psFileName_debut,"%s%s",psFileName,debut);
  e=sprintf(psFileName_fin,"%s%s",psFileName,fin);
  //	cout<<"debut "<<psFileName_debut<<endl;
  //	cout<<"fin "<<psFileName_fin<<endl;

  // Create a canvas into which we will draw our plots.
  Int_t canvas_w = 420;
  Int_t canvas_h = (Int_t)canvas_w*sqrt(2);
  TCanvas* c1 = new TCanvas("c1","n Cross-section",canvas_w,canvas_h);
  // Set logarithmic y axis
  c1->SetLogy();
  // Ticks up-right
  c1->SetTicks();
  c1->SetGrid();
  c1->cd();

  // Prepare the cuts and normalisation factors
  Double_t const pi=3.141592654;
  Double_t emin  = 1.;
  Double_t emax  = 1500.;
  Int_t    logE  = 1; // log E scale True=1, False=0
  Double_t fnorm = 3792.89/(100000.0*2.0*pi);

  cout << "Normalisation factor: " << fnorm << endl;
	
  TCut neutron = "(Avv==1)&&(Zvv==1)";

  TCut thet0 = "Tetlab<2.5";
  Double_t fnor0=fnorm/(1.-cos(2.5*pi/180.));
  TCut thet10 = "Tetlab>7.5 && Tetlab<12.5";
  Double_t fnor10=fnorm/(cos(7.5*pi/180.)-cos(12.5*pi/180.));
  TCut thet25 = "Tetlab>22.5 && Tetlab<27.5";
  Double_t fnor25=fnorm/(cos(22.5*pi/180.)-cos(27.5*pi/180.));
  TCut thet40 = "Tetlab>37.5 && Tetlab<42.5";
  Double_t fnor40=fnorm/(cos(37.5*pi/180.)-cos(42.5*pi/180.));
  TCut thet55 = "Tetlab>52.5 && Tetlab<57.5";
  Double_t fnor55=fnorm/(cos(52.5*pi/180.)-cos(57.5*pi/180.));
  TCut thet70 = "Tetlab>67.5 && Tetlab<72.5";
  Double_t fnor70=fnorm/(cos(67.5*pi/180.)-cos(72.5*pi/180.));
  TCut thet85 = "Tetlab>82.5 && Tetlab<87.5";
  Double_t fnor85=fnorm/(cos(82.5*pi/180.)-cos(87.5*pi/180.));
  TCut thet100 = "Tetlab>97.5 && Tetlab<102.5";
  Double_t fnor100=fnorm/(cos(97.5*pi/180.)-cos(102.5*pi/180.));
  TCut thet115 = "Tetlab>112.5 && Tetlab<117.5";
  Double_t fnor115=fnorm/(cos(112.5*pi/180.)-cos(117.5*pi/180.));
  TCut thet130 = "Tetlab>127.5 && Tetlab<132.5";
  Double_t fnor130=fnorm/(cos(127.5*pi/180.)-cos(132.5*pi/180.));
  TCut thet145 = "Tetlab>142.5 && Tetlab<147.5";
  Double_t fnor145=fnorm/(cos(142.5*pi/180.)-cos(147.5*pi/180.));
  TCut thet160 = "Tetlab>157.5 && Tetlab<162.5";
  Double_t fnor160=fnorm/(cos(157.5*pi/180.)-cos(162.5*pi/180.));

  // Plot theoretical cross sections:
  TH1F* hist0;
  Int_t first=1;
  Double_t fnora=1.;
  //plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet0,fnor0,fnora,"0 Deg",first);
  //first=0;
  first=1;
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet10,fnor10,fnora,"10^{o}",first);
  first=0;
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet10,fnor10,fnora,"10^{o}",first, kBlack);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet25,fnor25,fnora,"25^{o}",first);
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet25,fnor25,fnora,"25^{o}",first, kBlack);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet40,fnor40,fnora,"40^{o}",first);
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet40,fnor40,fnora,"40^{o}",first, kBlack);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet55,fnor55,fnora,"55^{o}",first);
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet55,fnor55,fnora,"55^{o}",first, kBlack);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet70,fnor70,fnora,"70^{o}",first);
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet70,fnor70,fnora,"70^{o}",first, kBlack);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet85,fnor85,fnora,"85^{o}",first);
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet85,fnor85,fnora,"85^{o}",first, kBlack);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet100,fnor100,fnora,"100^{o}",first);
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet100,fnor100,fnora,"100^{o}",first, kBlack);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet115,fnor115,fnora,"115^{o}",first);
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet115,fnor115,fnora,"115^{o}",first, kBlack);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet130,fnor130,fnora,"130^{o}",first);
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet130,fnor130,fnora,"130^{o}",first, kBlack);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet145,fnor145,fnora,"145^{o}",first);
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet145,fnor145,fnora,"145^{o}",first, kBlack);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c1,emin,emax,logE,neutron,thet160,fnor160,fnora,"160^{o}",first);
  plotTheorTheta(reffort,titre,c1,emin,emax,logE,neutron,thet160,fnor160,fnora,"160^{o}",first, kBlack);

  // Experimental points:
  Float_t fnorexp=1.;
  //  plotExpTheta("p1200pb20_000",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_010",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_025",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_040",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_055",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_070",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_085",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_100",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_115",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_130",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_145",racine,fnorexp);
//   fnorexp = fnorexp/10.;
//   plotExpTheta("p1200pb20_160",racine,fnorexp);


  // Legende
  //TLegend* legend = new TLegend(0.5,0.82,0.9,0.9, "");
  //  legend->SetFillColor(0);
  // legend->SetTextSize(0.0204499);
  //legend->AddEntry(hist0, "Reference INCL4+ABLA", "l");
  //legend->Draw("same");

  // The title (hist0->SetTitle(titre); also needed)!!!!
  c1->SetFillColor(0);
  TPaveText* pt = new TPaveText(0.153425,0.916155,0.854795,0.965235,"blNDC");
  pt->SetName("title");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetTextSize(0.04);
  pt->AddText(titre);
  pt->Draw();
   
  TTimeStamp* tt= new TTimeStamp(); // Date and time
  //	Char_t* date_time= tt->AsString("l");
  //	cout<<"time and date: "<<date_time<<endl;
  //   TPaveText* dt = new TPaveText(0.95,0.98,1.0,1.0,"NDC");
  //   dt->AddText(tt->AsString("l"));
  //   dt->SetBorderSize(0);
  //   dt->SetFillColor(0);
  //   dt->SetTextSize(0.015);
  //   dt->Draw();
   
  c1->Modified();
  c1->cd();
  c1->SetSelected(c1);
   
  c1->Print(psFileName_debut,"Portrait");


  // Second page of graphics:

  // Create a canvas into which we will draw our plots.
  //    Int_t canvas_w=420;
  //	Int_t canvas_h=canvas_w*sqrt(2);
  TCanvas* c2 = new TCanvas("c2","n Cross-section",canvas_w,canvas_h);

  c2->SetLogy();
  c2->SetTicks();
  c2->SetGrid();
  c2->cd();

  // Prepare the cuts and normalisation factors
  emin=0.;
  emax=1300.;
  logE=0; // log E scale True=1, False=0

  // Plot theoretical cross sections:
  first=1;
  fnora=1.;
  plotTheorTheta(ref,titre,c2,emin,emax,logE,neutron,thet0,fnor0,fnora,"0 Deg",first);
  first=0;
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c2,emin,emax,logE,neutron,thet10,fnor10,fnora,"10 Deg",first);
  fnora=fnora/10.;
  plotTheorTheta(ref,titre,c2,emin,emax,logE,neutron,thet25,fnor25,fnora,"25 Deg",first);

  // Experimental points:
  fnorexp=1.;
  //  plotExpTheta("p1200pb20_000",racine,fnorexp);
  fnorexp = fnorexp/10.;
  //  plotExpTheta("p1200pb20_010",racine,fnorexp);
  fnorexp = fnorexp/10.;
  //  plotExpTheta("p1200pb20_025",racine,fnorexp);

  //	legend->Draw("same");
  pt->Draw();
  //   	dt->Draw();
	
  c2->Modified();
  c2->cd();
  c2->SetSelected(c2);
   
  c2->Print(psFileName,"Portrait");

  // Third page of graphics:

  // Create a canvas into which we will draw our plots.
  //	Int_t canvas_w=420;
  //	Int_t canvas_h=canvas_w*sqrt(2);
  TCanvas* c3 = new TCanvas("c3","n Cross-section",canvas_w,canvas_h);
  c3->Range(0,0,1,1);
  TPad *c3_1 = new TPad("c3_1", "E*",0.,0.5,1.0,1.0);

  c3_1->SetLogy();
  c3_1->SetTicks();
  c3_1->SetGrid();
  c3_1->Draw();
  c3_1->cd();
  emin=0.;
  emax=1200.;
  Int_t nbx=100;
 
  TH1F* hista=new TH1F("hista","",nbx,emin,emax);
  hista->SetFillStyle(0);
  hista->SetStats(kTRUE);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  hista->SetTitle(titre);
  hista->GetXaxis()->SetTitle("Excitation Energy (MeV)");
  hista->GetXaxis()->CenterTitle(true);
  hista->GetXaxis()->SetLabelSize(0.03);
  hista->GetYaxis()->SetTitle("Counts");
  hista->GetYaxis()->CenterTitle(true);
  hista->GetYaxis()->SetLabelSize(0.03);
  hista->GetYaxis()->SetTitleOffset(1.3);
  hista->SetLineColor(kRed);
  ref->Draw("Exini >> hista");
  hista->DrawCopy(); // Needed to delete hista after

  //   dt->Draw();
  hista->Delete();

  c3->cd();
   
  TPad* c3_2 = new TPad("c3_2", "",0.,0.,1.0,0.5);
  c3_2->Draw();
  c3_2->cd();
  c3_2->Divide(2,1);
  TVirtualPad* c3_2_1 = c3_2->GetPad(1);   
  TVirtualPad* c3_2_2 = c3_2->GetPad(2);   

  c3_2_1->SetLogy();
  c3_2_1->SetTicks();
  c3_2_1->SetGrid();
  c3_2_1->cd();

  Double_t amax=210.5;
  Double_t amin=159.5;
  nbx=51;
  TH1F* hist1=new TH1F("hist1","",nbx,amin,amax);
  hist1->SetFillStyle(0);
  hist1->SetStats(kTRUE);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.9);
  hist1->SetTitle(titre);
  hist1->GetXaxis()->SetTitle("Remnant mass (A)");
  hist1->GetXaxis()->CenterTitle(true);
  hist1->GetXaxis()->SetLabelSize(0.03);
  hist1->GetYaxis()->SetTitle("Counts");
  hist1->GetYaxis()->CenterTitle(true);
  hist1->GetYaxis()->SetLabelSize(0.03);
  hist1->GetYaxis()->SetTitleOffset(1.3);
  hist1->SetLineColor(kRed);
  ref->Draw("Massini >> hist1");
  hist1->DrawCopy();
  hist1->Delete();

  c3_2->cd();
  c3_2_2->SetLogy();
  c3_2_2->SetTicks();
  c3_2_2->SetGrid();
  c3_2_2->cd();

  Double_t zmax=90.5;
  Double_t zmin=59.5;
  nbx=31;
  TH1F* hist2=new TH1F("hist2","",nbx,zmin,zmax);
  hist2->SetFillStyle(0);
  hist2->SetStats(kTRUE);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.9);
  hist2->SetTitle(titre);
  hist2->GetXaxis()->SetTitle("Remnant charge (Z)");
  hist2->GetXaxis()->CenterTitle(true);
  hist2->GetXaxis()->SetLabelSize(0.03);
  hist2->GetYaxis()->SetTitle("Counts");
  hist2->GetYaxis()->CenterTitle(true);
  hist2->GetYaxis()->SetLabelSize(0.03);
  hist2->GetYaxis()->SetTitleOffset(1.3);
  hist2->SetLineColor(kRed);
  ref->Draw("Mzini >> hist2");
  hist2->DrawCopy();
  hist2->Delete();

  //   dt->Draw();
  c3->Print(psFileName_fin,"Portrait");

  return;
}
