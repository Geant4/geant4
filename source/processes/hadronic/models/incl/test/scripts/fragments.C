// Procedure for residual nuclei cross sections   AB 10/2007
// Requires: test/tmp/run2.root datafile (generate by ./run.sh run2)
// Original version: Alain Boudard
// Modified by: Pekka Kaitaniemi 23.1.2008: code cleanup for Geant4 release
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "Riostream.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"
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

void fragments()
{
  gROOT->Reset();
  gROOT->ProcessLine(".x scripts/rootlogon.C");
  residus();
  // Neutron kinetic energy. Phys.Rev. article p21

  INCL(); // Plot INCL output variables: Remnant mass, charge, excitation and impact parameter
}

void INCL()
{
  TFile *cf = new TFile("tmp/run1.root");
  TFile *ff = new TFile("tmp/run1ref.root");

  TTree *ct = (TTree *) cf->Get("h101");
  ct->SetLineColor(kRed);
  TTree *ft = (TTree *) ff->Get("h101");

  TCanvas *inclResult = new TCanvas();
  inclResult->Divide(2,2);

  inclResult->cd(1);
  ft->Draw("Massini");
  ct->Draw("Massini", "", "same");

  inclResult->cd(2);
  ft->Draw("Mzini");
  ct->Draw("Mzini", "", "same");

  inclResult->cd(3);
  ft->Draw("Exini");
  ct->Draw("Exini", "", "same");

  inclResult->cd(4);
  ft->Draw("Jremn");
  ct->Draw("Jremn", "", "same");

  TCanvas *bimpactPlot = new TCanvas();
  ft->Draw("Bimpact");
  ct->Draw("Bimpact", "", "same");

  inclResult->SaveAs("run2INCL.eps");
}

void residus() {
  TFile* fref = new TFile("tmp/run1.root");
  TFile* freffort = new TFile("tmp/run1ref.root");
  Double_t sigGeom1 = 3792.89;
  Double_t nbEvt1 = 100000.0;
	
  Char_t* titre="p(1.0 GeV) + 208Pb (INCL4+ABLA)";
	
  // Graphic extensions (Z and A):
  Float_t  zemin=5;
  Float_t  zemax=85;
  Double_t sigZmin=1.e+0;
  Double_t sigZmax=1.e+3; 
  Float_t  aemin=10;
  Float_t  aemax=210;
  Double_t sigAmin=2.e-2;
  Double_t sigAmax=2.e+2;
		
  Char_t* psFileName="GSI_residues.ps"; // PS file produced
	
  // The path to experimental files and experimental file name:
  Char_t* racine="./data/data_gsi/pb_p_1000/";
  Char_t* file_exp="spall_pb.txt";
  Char_t* file_exp2="fiss_pb.txt";

  // Prepare the normalisation factors
  Double_t fnorm1 = sigGeom1/nbEvt1;
  cout << "Normalisation factor: "<<fnorm1 << endl;
  
  TTree* ref = (TTree*) fref->Get("h101"); // Read the ROOT tree
  TTree* reffort = (TTree*) freffort->Get("h101");
  
  // Prepare the postscript file name for several pages:
  Char_t psFileName_debut[256];
  Char_t psFileName_fin[256];
  Int_t e;
  Char_t* debut ="(";
  Char_t* fin =")";
  e=sprintf(psFileName_debut,"%s%s",psFileName,debut);
  e=sprintf(psFileName_fin,"%s%s",psFileName,fin);

  // Experimental points:
  Char_t file_exp_path[256];
  Char_t file_exp_path2[256];

  //Concatenation Racine+File_Name
  e=sprintf(file_exp_path,"%s%s",racine,file_exp);
  cout<<"Experimental file evaporation: "<<file_exp_path<<endl;
  e=sprintf(file_exp_path2,"%s%s",racine,file_exp2);
  cout<<"Experimental file fission: "<<file_exp_path2<<endl;
	
  // Read experimental data (2 files for Lead, fission and evapo):
  FILE* fexp =  fopen(file_exp_path,"r");
  const Int_t dimexp=1000;
  Float_t zz[dimexp];
  Float_t aa[dimexp];
  Float_t sig[dimexp];
  Float_t dsig[dimexp];
	
  Int_t iexp=0;
  do {
    e=fscanf(fexp,"%f %f %f %f ",&zz[iexp],&aa[iexp],&sig[iexp],&dsig[iexp]);
    if(e==EOF) break;
    //	cout<<zz[iexp]<<"  "<<aa[iexp]<<"  "<<sig[iexp]<<"  "<<dsig[iexp]<<endl;
    iexp=iexp+1;
    if(iexp>dimexp) {cout<<"More experimental values than the dimension"<<endl;
    return;}
  }
  while(1);
  FILE* fexp2 = fopen(file_exp_path2,"r");
  do {
    e=fscanf(fexp2,"%f %f %f %f ",&zz[iexp],&aa[iexp],&sig[iexp],&dsig[iexp]);
    if(e==EOF) break;
    //	cout<<zz[iexp]<<"  "<<aa[iexp]<<"  "<<sig[iexp]<<"  "<<dsig[iexp]<<endl;
    iexp=iexp+1;
    if(iexp>dimexp) {cout<<"More experimental values than the dimension"<<endl;
    return;}
  }
  while(1); 

  // Z and A cross sections
  const Int_t dimZ=100;
  const Int_t dimA=300;
  Float_t sigZExp[dimZ];
  Float_t dsigZExp[dimZ];
  Float_t zexp[dimZ];
  Float_t dzexp[dimZ];
  Float_t sigAExp[dimA];
  Float_t dsigAExp[dimA];
  Float_t aexp[dimA];
  Float_t daexp[dimA];
  Int_t Z;
  Int_t A;

  // INIT and Put tables to 0 for possible multi calls in root.
  for(Int_t i=0;i<dimZ;++i) {
    dsigZExp[i] = 0.;
    sigZExp[i]  = 0.;
    zexp[i]     = (Float_t)i;
    dzexp[i]    = 0.;
  }
  for(Int_t i=0;i<dimA;++i) {
    dsigAExp[i] = 0.;
    sigAExp[i]  = 0.;
    aexp[i]     = (Float_t)i;
    daexp[i]    = 0.;
  }
  TH2D* verite= new TH2D("verite","",dimZ,-0.5,99.5,dimA,-0.5,299.5);

  for(Int_t i=0;i<iexp;++i) { // Sum cross sections
    Z=(Int_t)zz[i];
    A=(Int_t)aa[i];
    sigZExp[Z]  = sigZExp[Z]  + sig[i];
    dsigZExp[Z] = dsigZExp[Z] + dsig[i]*dsig[i];
    sigAExp[A]  = sigAExp[A]  + sig[i];
    dsigAExp[A] = dsigAExp[A] + dsig[i]*dsig[i];

    // Here Bin numbers and not Z,A values (Z+1 due to the choice dimZ,-0.5,99.5)
    verite->SetBinContent(Z+1,A+1,fnorm1);
  }

  // Experimental errors of summed cross sections:	
  for(Int_t i=0;i<dimZ;++i) {dsigZExp[i]=sqrt(dsigZExp[i]);}
  for(Int_t i=0;i<dimA;++i) {dsigAExp[i]=sqrt(dsigAExp[i]);}
	
  // Create a canvas into which we will draw our plots.
  Int_t canvas_w=420;
  Int_t canvas_h=(Int_t)canvas_w*sqrt(2);
  TCanvas* c1 = new TCanvas("c1","A and Z residue cross-sections",canvas_w,canvas_h);
  c1->Range(0,0,1,1);
  TPad *c1_1 = new TPad("c1_1", "Z",0.,0.5,1.0,1.0);
  TPad *c1_2 = new TPad("c1_2", "A",0.,0.,1.0,0.5);
	
  // Sum theory on measured elementary cross sections and normalize:	
  TH2D* calc1= new TH2D("calc1","",dimZ,-0.5,99.5,dimA,-0.5,299.5);
  ref->Draw("Avv:Zvv >> calc1");
  calc1->Multiply(calc1,verite);

  TH2D* calc1_fort= new TH2D("calc1_fort","",dimZ,-0.5,99.5,dimA,-0.5,299.5);
  reffort->Draw("Avv:Zvv >> calc1_fort");
  calc1_fort->Multiply(calc1_fort,verite);
	
  // Z Plot
  c1_1->SetLogy();
  c1_1->SetTicks();
  c1_1->SetGrid();
  c1_1->Draw();
  c1_1->cd();
  //	Float_t zmin=-0.5;
  Float_t zmin=zemin-0.5;
  Float_t zmax=zemax+0.5;
  Int_t nbx=(Int_t)(zmax-zmin);
  //	cout<<zmax<<"  "<<zmin<<"  "<<nbx<<endl;
	
  TH1D* histz_on_measured=new TH1D("histz_on_measured","",nbx,zmin,zmax);
  histz_on_measured=calc1->ProjectionX();
	
  //	histz_on_measured->Print("all");
  //	for(Int_t i=1;i<nbx+1;++i) cout<<i<<"  "<<histz_on_measured->GetBinContent(i)<<" projected (bin number)"<<endl;
	
  histz_on_measured->SetLineColor(kRed);
 
  TH1F* histz=new TH1F("histz","",nbx,zmin,zmax);
  //	cout<<"histz "<<nbx<<"  "<<zmin<<"  "<<zmax<<endl;
  histz->SetFillStyle(0);
  histz->SetStats(kFALSE);
  histz->SetTitle(titre);
  histz->GetXaxis()->SetTitle("Z of the residue");
  histz->GetXaxis()->CenterTitle(true);
  histz->GetXaxis()->SetLabelSize(0.03);
  histz->GetYaxis()->SetTitle("Cross-Section (mb)");
  histz->GetYaxis()->CenterTitle(true);
  histz->GetYaxis()->SetLabelSize(0.03);
  histz->GetYaxis()->SetTitleOffset(1.3);
  histz->SetLineColor(kGreen);
		
  ref->Draw("Zvv >> histz");
  histz->Scale(fnorm1);
	
  // Full y extension of the picture:
  histz->SetMinimum(sigZmin);
  histz->SetMaximum(sigZmax);
	
  histz->DrawCopy();	// Needed to delete histz after
  histz_on_measured->DrawCopy("Same");
	
  TGraphErrors* grexpz= new TGraphErrors(dimZ,zexp,sigZExp,dzexp,dsigZExp);
  grexpz->SetMarkerColor(4);
  grexpz->SetMarkerStyle(21);
  grexpz->SetMarkerSize(0.4);

  grexpz->Draw("PZ");
	
  TTimeStamp* tt= new TTimeStamp(); // Date and time
  TPaveText* dt = new TPaveText(0.95,0.96,1.0,1.0,"NDC");
  dt->AddText(tt->AsString("l"));
  dt->SetBorderSize(0);
  dt->SetFillColor(0);
  dt->SetTextSize(0.03);
  dt->Draw();

  // A plot
  c1->cd();
  c1_2->SetLogy();
  c1_2->SetTicks();
  c1_2->SetGrid();
  c1_2->Draw();
  c1_2->cd();
  Float_t amin=aemin-0.5;
  Float_t amax=aemax+0.5;
  nbx=(Int_t)(amax-amin);
	
  TH1D* hista_on_measured=new TH1D("hista_on_measured","",nbx,amin,amax);
  hista_on_measured=calc1->ProjectionY();
  hista_on_measured->SetLineColor(kRed);

  TH1D* hista_on_measured_fortran=new TH1D("hista_on_measured_fortran","",nbx,amin,amax);
  hista_on_measured_fortran=calc1_fort->ProjectionY();
  hista_on_measured_fortran->SetLineColor(kBlack);
 
  TH1F* hista=new TH1F("hista","",nbx,amin,amax);
  hista->SetFillStyle(0);
  hista->SetStats(kFALSE);
  hista->GetXaxis()->SetTitle("A of the residue");
  hista->GetXaxis()->CenterTitle(true);
  hista->GetXaxis()->SetLabelSize(0.03);
  hista->GetYaxis()->SetTitle("Cross-Section (mb)");
  hista->GetYaxis()->CenterTitle(true);
  hista->GetYaxis()->SetLabelSize(0.03);
  hista->GetYaxis()->SetTitleOffset(1.3);
  hista->SetLineColor(kRed);
		
  ref->Draw("Avv >> hista");
  hista->Scale(fnorm1);
	
  // Full y extension of the picture:
  hista->SetMinimum(sigAmin);
  hista->SetMaximum(sigAmax);
	
  //hista->DrawCopy();	// Needed to delete hista after
  hista_on_measured->DrawCopy("Same");
  //hista_on_measured->DrawCopy();
  TH1F* histfortran=new TH1F("histfortran","",nbx,amin,amax);
  histfortran->SetLineColor(kBlack);
  reffort->SetLineColor(kBlack);
  reffort->Draw("Avv >> histfortran", "", "same");
  histfortran->Scale(fnorm1);

  histfortran->Draw("same");

  TGraphErrors* grexpa= new TGraphErrors(dimA,aexp,sigAExp,daexp,dsigAExp);
  grexpa->SetMarkerColor(4);
  grexpa->SetMarkerStyle(21);
  grexpa->SetMarkerSize(0.4);

  grexpa->Draw("PZ");
	 
  //c1->Print(psFileName_debut,"Portrait"); // saving the ps file

  c1->Clear();
  c1->SetLogy(1);
  hista->SetTitle(titre);
  hista_on_measured->Draw();
  hista->Draw("same");
  histfortran->Draw("same");
  grexpa->Draw("PZ");
  //  c1->SaveAs("fragments.png");
  //  c1->SaveAs("fragments.eps");
  // c1->Print(psFileName,"Portrait"); // saving the ps file

  c1->Clear();
  c1->SetLogy(1);
  hista_on_measured->GetXaxis()->SetTitle("Mass number");
  hista_on_measured->GetYaxis()->SetTitle("Cross section (mb)");
  hista_on_measured->SetTitle(titre);
  hista_on_measured->Draw();
  hista_on_measured_fortran->Draw("same");
  grexpa->Draw("PZ");
  c1->SaveAs("fragments.png");
  c1->SaveAs("fragments.eps");
  c1->Print(psFileName_debut,"Portrait"); // saving the ps file

  TCanvas *remnantCanvas = new TCanvas();
  remnantCanvas->Divide(2,1);
  remnantCanvas->cd(1);
  remnantCanvas->SetLogy(1);
  ref->SetLineColor(kRed);
  reffort->Draw("Massini");
  ref->Draw("Massini", "", "same");
  reffort->GetHistogram()->SetTitle("p(1 GeV) + 208Pb after INCL cascade");
  reffort->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
  reffort->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
  reffort->GetHistogram()->GetXaxis()->SetTitle("Remnant nucleus mass number");
  reffort->GetHistogram()->GetYaxis()->SetTitle("Counts");
  remnantCanvas->cd(2);
  remnantCanvas->SetLogy(1);
  reffort->Draw("Exini", "", "");
  ref->Draw("Exini", "", "same");
  reffort->GetHistogram()->SetTitle("");
  reffort->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
  reffort->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
  reffort->GetHistogram()->GetXaxis()->SetTitle("Remnant nucleus excitation energy");
  reffort->GetHistogram()->GetYaxis()->SetTitle("Counts");
  remnantCanvas->SaveAs("run2remnants.png");

  // Second page, isotopic cross sections
  Int_t zDebut=84;
  Int_t zFin=20;
  Int_t zCourant=zDebut;
	
  // Name of nucleus:
  Char_t* NucName="   ";
  TPaveLabel* pNucName;
  extern Char_t* fNucName(Int_t);
 
  // Variables for experimental cross sections:
  const Int_t nbExpIsoXS = 30;
  Float_t sigIso[nbExpIsoXS];
  Float_t dsigIso[nbExpIsoXS];
  Float_t AIso[nbExpIsoXS];
  Float_t BidA[nbExpIsoXS];
  Int_t j;
  TGraphErrors* griso;
	
  // Several pages of residues:
  Int_t flag_NoPlot=0;
  TCanvas* ciso;
  Char_t canvaName[128];
  Int_t l=0;
  do {
    // saving the ps file:   
    if (l!=0) ciso->Print(psFileName,"Portrait");
	
    l=l+1;	//page number
    cout<<"l= "<<l<<endl;
    switch (l) {
    case 1:
      amin=149.5;	//extensions for this page
      amax=210.5;
      sigAmin=2.e-2;
      sigAmax=2.e+2;
      break;
    case 2:
      amin=109.5;
      amax=170.5;
      sigAmin=1.e-3;
      sigAmax=1.e+1;
      break;
    case 3:
      amin=69.5;
      amax=130.5;
      sigAmin=5.e-3;
      sigAmax=0.5e+1;
      break;
    case 4:
      amin=39.5;
      amax=100.5;
      sigAmin=1.e-2;
      sigAmax=1.e+1;
      break;
    case 5:
      amin=4.5;
      amax=65.5;
      sigAmin=1.e-2;
      sigAmax=1.e+1;
      break;
    default:
      amin=149.5;	//default extensions
      amax=210.5;
      sigAmin=2.e-2;
      sigAmax=2.e+2;
    }

    // Create a canvas into which we will draw our plots.
    canvas_w=420;
    canvas_h=(Int_t)canvas_w*sqrt(2);
    e=sprintf(canvaName,"ciso%d",l);
    cout<<"Canva name "<<canvaName<<endl;
    ciso = new TCanvas(canvaName,"Isotopic cross sections",canvas_w,canvas_h);
    ciso->Range(0,0,1,1);
    ciso->cd();
	
    // Prepare subpad dimensions:	
    Float_t margeX=0.05;
    Float_t margeY=0.1;
    const Int_t nbGraphesX=3;
    const Int_t nbGraphesY=5;
	
    const Int_t nbGraphes=nbGraphesX*nbGraphesY;
    const Int_t nbGraphesMax=30;
	
    if(nbGraphes>nbGraphesMax) { 
      cout<<"Number of subGraphes must be smaller than nbGraphesMax= "<<nbGraphesMax<<endl;
      return;
    }
	
    Float_t unitWidthX=(1.-margeY)/nbGraphesX;
    Float_t unitWidthY=(1.-margeX)/nbGraphesY;
    TPad* c2[nbGraphes];
    Char_t c2name[nbGraphes];
    Char_t* c2char="c2_";
	
    nbx=(Int_t)(amax-amin);
    TH1F* histiso=new TH1F("histiso","",nbx,amin,amax);
    histiso->SetFillStyle(0);
    histiso->SetStats(kFALSE);
    histiso->GetXaxis()->SetTitle("A (residues)");
    histiso->GetXaxis()->SetTitleSize(0.1);
    histiso->GetXaxis()->CenterTitle(true);
    histiso->GetXaxis()->SetLabelSize(0.08);
    histiso->GetYaxis()->SetTitle(" #sigma (mb)");
    histiso->GetYaxis()->SetTitleSize(0.1);
    histiso->GetYaxis()->CenterTitle(true);
    histiso->GetYaxis()->SetLabelSize(0.1);
    histiso->GetYaxis()->SetTitleOffset(1.0);
    histiso->SetLineColor(kRed);

    TH1F* histisofort=new TH1F("histisofort","",nbx,amin,amax);
    histisofort->SetFillStyle(0);
    histisofort->SetStats(kFALSE);
    histisofort->GetXaxis()->SetTitle("A (residues)");
    histisofort->GetXaxis()->SetTitleSize(0.1);
    histisofort->GetXaxis()->CenterTitle(true);
    histisofort->GetXaxis()->SetLabelSize(0.08);
    histisofort->GetYaxis()->SetTitle(" #sigma (mb)");
    histisofort->GetYaxis()->SetTitleSize(0.1);
    histisofort->GetYaxis()->CenterTitle(true);
    histisofort->GetYaxis()->SetLabelSize(0.1);
    histisofort->GetYaxis()->SetTitleOffset(1.0);
    histisofort->SetLineColor(kBlack);

    Int_t ix; Int_t iy;
    Float_t X1; Float_t Y1; Float_t X2; Float_t Y2;
    Char_t zSelect[nbGraphes];

    for(Int_t i=0;i<=nbGraphes-1;++i) {
      ix=i%3;
      iy=i/3;
      e=sprintf(zSelect,"Zvv==%d",zCourant);
      cout<<"ZCourant= "<<zCourant<<" select= "<<zSelect<<endl;
      e=sprintf(c2name,"%s%d",c2char,i);

      if(ix==0) {
	X1 = ix*unitWidthX;
	X2 = X1+unitWidthX + margeY; //first graph of the line
      }	else {
	X1 = ix*unitWidthX + margeY;
	X2 = X1+unitWidthX;
      }

      if(iy==(nbGraphesY-1)) {
	Y1 = (nbGraphesY-iy-1)*unitWidthY;
	Y2 = Y1+unitWidthY + margeX; //last graph of the column
      } else {
	Y1 = (nbGraphesY-iy-1)*unitWidthY + margeX;
	Y2 = Y1+unitWidthY;
      }
		
      c2[i]=new TPad(c2name,"",X1,Y1,X2,Y2);
		
      if(ix==0) {
	c2[i]->SetLeftMargin(margeY/(unitWidthX + margeY));
      } else {
	c2[i]->SetLeftMargin(0.01);
      }

      if(iy==(nbGraphesY-1)) {
	c2[i]->SetBottomMargin(margeX/(unitWidthY + margeX));
      } else {
	c2[i]->SetBottomMargin(0.01);
      }

      c2[i]->SetFillColor(0);
      c2[i]->SetBorderMode(0);
      c2[i]->SetBorderSize(2);
      c2[i]->SetRightMargin(0.01);
      c2[i]->SetTopMargin(0.01);
      c2[i]->SetFrameBorderMode(0);
      c2[i]->SetFrameBorderMode(0);
      c2[i]->SetLogy();
      c2[i]->SetTicks();
      c2[i]->SetGrid();
      c2[i]->Draw();
      c2[i]->cd();
	
      if(flag_NoPlot==0) {
	ref->Draw("Avv >> histiso",zSelect);
	reffort->Draw("Avv >> histisofort",zSelect, "same");
      }
      else histiso->Reset();
	
      histiso->Scale(fnorm1);
      histisofort->Scale(fnorm1);
	
      histiso->SetMinimum(sigAmin);
      histiso->SetMaximum(sigAmax);
      histisofort->SetMinimum(sigAmin);
      histisofort->SetMaximum(sigAmax);
		
      histisofort->DrawCopy();
      histiso->DrawCopy("same");	// Needed to delete histiso after

      NucName = fNucName(zCourant);
      pNucName = new TPaveLabel(amin,sigAmax/4.,amin+(amax-amin)/3.,sigAmax,NucName);
      pNucName->SetBorderSize(0);
      pNucName->SetFillColor(0);
      pNucName->Draw("Same");
      //		ciso->cd();

      if(flag_NoPlot==0) {
	// Experimental values:
        j=0;
	for(Int_t k=0;k<iexp;++k) {
	  Z=(Int_t)zz[k];
	  //		cout<<Z<<" "<<zCourant<<endl;
	  if (Z!=zCourant) continue;
	  A=(Int_t)aa[k];
	  sigIso[j]=sig[k];
	  dsigIso[j]=dsig[k];
	  AIso[j]=aa[k];
	  BidA[j]=0.;
	  j=j+1;
	  if(j>nbExpIsoXS) {cout<<"enlarge nbExpIsoXS! "<<nbExpIsoXS<<endl; return;}
	}
	if (j!=0) {griso= new TGraphErrors(j,AIso,sigIso,BidA,dsigIso);
	griso->SetMarkerColor(4);
	griso->SetMarkerStyle(21);
	griso->SetMarkerSize(0.4); 
	griso->Draw("PZ");}
      }	   
		   
      ciso->cd();
      zCourant=zCourant-1;
      if(zCourant<zFin) {flag_NoPlot=1;}
    }
    ciso->cd();}	   
	
  while (zCourant>zFin);
	   
  ciso->Print(psFileName_fin,"Portrait"); // saving the ps file
	
  return;
}

Char_t* fNucName(Int_t Z) {
  Char_t* NucName[101];
  NucName[0]="n";
  if (Z>100) {cout<<" Enlarge NucName table!"<<endl; return NucName[0];}
  NucName[1]="p";
  NucName[2]="_{2}He";
  NucName[3]="_{3}Li";
  NucName[4]="_{4}Be";
  NucName[5]="_{5}B";
  NucName[6]="_{6}C";
  NucName[7]="_{7}N";
  NucName[8]="_{8}O";
  NucName[9]="_{9}F";
  NucName[10]="_{10}Ne";
  NucName[11]="_{11}Na";
  NucName[12]="_{12}Mg";
  NucName[13]="_{13}Al";
  NucName[14]="_{14}Si";
  NucName[15]="_{15}P";
  NucName[16]="_{16}S";
  NucName[17]="_{17}Cl";
  NucName[18]="_{18}Ar";
  NucName[19]="_{19}K";
  NucName[20]="_{20}Ca";
  NucName[21]="_{21}Sc";
  NucName[22]="_{22}Ti";
  NucName[23]="_{23}V";
  NucName[24]="_{24}Cr";
  NucName[25]="_{25}Mn";
  NucName[26]="_{26}Fe";
  NucName[27]="_{27}Co";
  NucName[28]="_{28}Ni";
  NucName[29]="_{29}Cu";
  NucName[30]="_{30}Zn";
  NucName[31]="_{31}Ga";
  NucName[32]="_{32}Ge";
  NucName[33]="_{33}As";
  NucName[34]="_{34}Se";
  NucName[35]="_{35}Br";
  NucName[36]="_{36}Kr";
  NucName[37]="_{37}Rb";
  NucName[38]="_{38}Sr";
  NucName[39]="_{39}Y";
  NucName[40]="_{40}Zr";
  NucName[41]="_{41}Nb";
  NucName[42]="_{42}Mo";
  NucName[43]="_{43}Tc";
  NucName[44]="_{44}Ru";
  NucName[45]="_{45}Rh";
  NucName[46]="_{46}Pd";
  NucName[47]="_{47}Ag";
  NucName[48]="_{48}Cd";
  NucName[49]="_{49}In";
  NucName[50]="_{50}Sn";
  NucName[51]="_{51}Sb";
  NucName[52]="_{52}Te";
  NucName[53]="_{53}I";
  NucName[54]="_{54}Xe";
  NucName[55]="_{55}Cs";
  NucName[56]="_{56}Ba";
  NucName[57]="_{57}La";
  NucName[58]="_{58}Ce";
  NucName[59]="_{59}Pr";
  NucName[60]="_{60}Nd";
  NucName[61]="_{61}Pm";
  NucName[62]="_{62}Sm";
  NucName[63]="_{63}Eu";
  NucName[64]="_{64}Gd";
  NucName[65]="_{65}Tb";
  NucName[66]="_{66}Dy";
  NucName[67]="_{67}Ho";
  NucName[68]="_{68}Er";
  NucName[69]="_{69}Tm";
  NucName[70]="_{70}Yb";
  NucName[71]="_{71}Lu";
  NucName[72]="_{72}Hf";
  NucName[73]="_{73}Ta";
  NucName[74]="_{74}W";
  NucName[75]="_{75}Re";
  NucName[76]="_{76}Os";
  NucName[77]="_{77}Ir";
  NucName[78]="_{78}Pt";
  NucName[79]="_{79}Au";
  NucName[80]="_{80}Hg";
  NucName[81]="_{81}Tl";
  NucName[82]="_{82}Pb";
  NucName[83]="_{83}Bi";
  NucName[84]="_{84}Po";
  NucName[85]="_{85}At";
  NucName[86]="_{86}Rn";
  NucName[87]="_{87}Fr";
  NucName[88]="_{88}Ra";
  NucName[89]="_{89}Ac";
  NucName[90]="_{90}Th";
  NucName[91]="_{91}Pa";
  NucName[92]="_{92}U";
  NucName[93]="_{93}Np";
  NucName[94]="_{94}Pu";
  NucName[95]="_{95}Am";
  NucName[96]="_{96}Cm";
  NucName[97]="_{97}Bk";
  NucName[98]="_{98}Cf";
  NucName[99]="_{99}Es";
  NucName[100]="_{100}Fm";
	
  return NucName[Z];
} 
