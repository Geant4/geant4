#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <list>

#include <math.h>
#include <vector>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TRefArray.h"
#include "TStyle.h"
#include "TGraph.h"

const int modelsITEP=9, modelsBNL=9, modelsD=9;

std::string ModelsITEP[9]  = {"lepar", "ftfb",    "bertini", "binary", "qgsc",      "qgsp",      "qgsb",   "ftfp",      "CHIPS"};
std::string ModelsITEPh[9] = {"lepar", "ftfb",    "bertini", "binary", "qgsc",      "qgsp",      "qgsb",   "ftfp",      "CHIPS"};
std::string ModelNamesI[9] = {"LEP",   "FTF-Bin", "Bertini", "Binary", "QGS-Chips", "QGS-Preco", "QGS-Bin","FTF-Preco", "CHIPS"};

std::string ModelNameD = "bertini";
std::string ModelDirectory[9] = {"9.3.cand05", "9.3.ref03", "V09-03-23", "9.3.p01", "V09-03-28", "9.3.ref05", "9.3.ref06", "V09-03-45", "V09-03-66"};
std::string ModelNamesD[9] = {"(9.3.cand05)", "(9.3.ref03)", "(V09-03-23)", "(9.3.p01)", "(V09-03-28)", "(9.3.ref05)", "(9.3.ref06)", "(V09-03-45)", "(V09-03-66)"};

std::string ModelsBNL[9]   = {"lepar", "ftfb",    "bertini", "binary", "qgsc",      "qgsp",      "qgsb",   "ftfp", "CHIPS"};
std::string ModelsBNLh[9]  = {"lepar", "ftfb",    "bertini", "binary", "qgsc",      "qgsp",      "qgsb",   "ftfp", "CHIPS"};
std::string ModelNamesB[9] = {"LEP",   "FTF-Bin", "Bertini", "Binary", "QGS-Chips", "QGS-Preco", "QGS-Bin","FTF-Preco","CHIPS"};


int         colModel[9]    = {8, 2, 6, 3, 7, 9, 1, 4, 12};
int         symbModel[9]   = {24, 29, 25, 27, 26, 23, 21, 20, 22};
int         stylModel[9]   = {3, 1, 2, 4, 5, 6, 1, 7, 8};
double      keproton[4]    = {0.09, 0.15, 0.19, 0.23};
double      keneutron[4]   = {0.07, 0.11, 0.15, 0.17};
bool        debug=false;

void Plot(int index, char test[5]="ITEP") {

  if ( test == "BNL" ) {
    switch(index) {
    case 1:
      plotMT4("Be","14.6", 0,1,1,-1.,-1., "piplus",  "proton");
      plotMT4("Be","14.6", 0,1,1,-1.,-1., "piminus", "proton");
      break;
    case 2:
      plotMT4("Cu","14.6", 0,1,1,-1.,-1., "proton", "proton");
      plotMT4("Cu","14.6", 0,1,1,-1.,-1., "kplus",  "proton");
      plotMT4("Cu","14.6", 0,1,1,-1.,-1., "kminus", "proton");
      break;
    case 3:
      plotMT4("Au","14.6", 0,1,1,-1.,-1., "piplus",  "proton");
      plotMT4("Au","14.6", 0,1,1,-1.,-1., "piminus", "proton");
      break;
    }
  } else if ( test == "ITEP" ) {
    TCanvas* cnv = new TCanvas("cnv","",400,300);
    cnv->cd(); cnv->SetLogy(1); cnv->SetLeftMargin(0.15);
    switch(index) {
    case 1:
      plotKE4("C","1.40",0,1,1,-1.,-1.,"proton","piplus");
      cnv->cd();
      plotKE("C","1.40","119.0",0,1,-1.,-1.,"neutron", "piplus");
      cnv->SaveAs("piplusCtoneutronat1.40GeV_1.eps");
      break;
    case 2:
      plotKE4("U","1.40",0,1,1,-1.,-1.,"proton","piplus");
      cnv->cd();
      plotKE("U","1.40","119.0",0,1,-1.,-1.,"neutron","piplus");
      cnv->SaveAs("piplusUtoneutronat1.40GeV_1.eps");
      break;
    case 3:
      plotKE4("C","5.00",0,1,1,-1.,-1.,"proton","piplus");
      cnv->cd();
      plotKE("C","5.00","119.0",0,1,-1.,-1.,"neutron","piplus");
      cnv->SaveAs("piplusCtoneutronat5.00GeV_1.eps");
      break;
    case 4:
      plotKE4("U","5.00",0,1,1,-1.,-1.,"proton","piplus");
      cnv->cd();
      plotKE("U","5.00","119.0",0,1,-1.,-1.,"neutron","piplus");
      cnv->SaveAs("piplusUtoneutronat5.00GeV_1.eps");
      break;
    case 5:
      plotKE4("C","5.00",0,1,1,-1.,-1.,"proton","piminus");
      cnv->cd();
      plotKE("C","5.00","119.0",0,1,-1.,"neutron","piminus");
      cnv->SaveAs("piminusCtoneutronat5.00GeV_1.eps");
      break;
    case 6:
      plotKE4("Cu","5.00",0,1,1,-1.,-1.,"proton","piminus");
      cnv->cd();
      plotKE("Cu","5.00","119.0",0,1,-1.,"neutron","piminus");
      cnv->SaveAs("piminusCutoneutronat5.00GeV_1.eps");
      break;
    case 7:
      plotKE4("Pb","5.00",0,1,1,-1.,-1.,"proton","piminus");
      cnv->cd();
      plotKE("Pb","5.00","119.0",0,1,-1.,-1.,"neutron","piminus");
      cnv->SaveAs("piminusPbtoneutronat5.00GeV_1.eps");
      break;
    case 8:
      plotKE4("U","5.00",0,1,1,-1.,-1.,"proton","piminus");
      cnv->cd();
      plotKE("C","5.00","119.0",0,1,-1.,-1.,"neutron","piminus");
      cnv->SaveAs("piminusUtoneutronat5.00GeV_1.eps");
      break;
    case 9:
      plotKE4("C","1.40",0,1,1);
      cnv->cd();
      plotKE("C","1.40","119.0",0,1,-1.,-1.,"neutron");
      cnv->SaveAs("protonCtoneutronat1.40GeV_1.eps");
      break;
    case 10:
      plotKE4("U","1.40",0,1,1);
      cnv->cd();
      plotKE("U","1.40","119.0",0,1,-1.,-1.,"neutron");
      cnv->SaveAs("protonUtoneutronat1.40GeV_1.eps");
      break;
    case 11:
      plotKE4("C","7.50",0,1,1);
      cnv->cd();
      plotKE("C","7.50","119.0",0,1,-1.,-1.,"neutron");
      cnv->SaveAs("protonCtoneutronat7.50GeV_1.eps");
      break;
    case 12:
      plotKE4("U","7.50",0,1,1);
      cnv->cd();
      plotKE("U","7.50","119.0",0,1,-1.,-1.,"neutron");
      cnv->SaveAs("protonUtoneutronat7.50GeV_1.eps");
      break;
    }
  }
  
  return;
}


void plotStandard(bool ratio=false, char dir[20]=".", char dird[40]=".", 
		  int leg1=1, int leg2=1, char mark=' ', bool error=true) {

  plotKEp("C", "5.00",0,1,1,-1.0,-1.,"piplus",ratio,error,leg1,leg2,dir,dird,mark);
  plotKEp("U", "5.00",0,1,1, 10.,-1.,"piplus",ratio,error,leg1,leg2,dir,dird,mark);
  plotKEp("C", "7.50",0,1,1, 1.0,-1.,"proton",ratio,error,leg1,leg2,dir,dird,mark);
  plotKEp("U", "7.50",0,1,1, 10.,-1.,"proton",ratio,error,leg1,leg2,dir,dird,mark);
  plotKEx("5.00", " 59.1", 0,1,1, -1.,-1., "proton",  "piminus", ratio,error,
	  leg1,leg2, dir,dird, mark);
  plotKEx("5.00", "119.0", 0,1,1,100.,-1., "neutron", "piminus", ratio,error,
	  leg1,leg2, dir,dird, mark);
  plotKEn("5.00", 0,1,1, 10.,-1., "piplus", ratio,error,leg1,leg2, dir,dird, mark);
  plotKEn("7.50", 0,1,1, 10.,-1., "proton", ratio,error,leg1,leg2, dir,dird, mark);
  plotMT4(     "14.6", 0,1,1,-1.,-1.,"piplus", "proton", ratio,error,
	  leg1,leg2, dir,dird, mark);
  plotMT4(     "14.6", 0,1,1,-1.,-1.,"piminus","proton",ratio,error,
          leg1,leg2, dir,dird, mark);
  plotMT4("Cu","14.6", 0,1,1,-1.,-1.,"kplus",  "proton",ratio,error,
	  leg1,leg2, dir,dird, mark);
  plotMT4("Cu","14.6", 0,1,1,-1.,-1.,"kminus", "proton",ratio,error,
	  leg1,leg2, dir,dird, mark);
  plotMT4("Cu","14.6", 0,1,1,-1.,-1.,"proton", "proton",ratio,error,
	  leg1,leg2, dir,dird, mark);
}


void plotKEx(char ene[6], char angle[6], int first=0, int logy=0, int save=0, 
	     double ymin=-1., double ymax=-1.,char particle[8]="proton", 
	     char beam[8]="proton", bool ratio='false', bool error=true, 
	     int leg1=1, int leg2=1, char dir[20]=".", char dird[40]=".", 
	     char mark=' ') {

  setStyle();
  TCanvas *myc = new TCanvas("myc","",800,600); myc->Divide(2,2);

  int leg=leg1; if (leg2 < leg) leg=leg2;
  char markf[4]=" ";
  myc->cd(1); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(a)");
  if (ratio)
    plotKERatio("C", ene,angle,first,logy,ymin,ymax,particle,beam,error,leg1,dir,dird,markf);
  else
    plotKE("C", ene,angle,first,logy,ymin,ymax,particle,beam, leg1, dir,dird,markf);
  myc->cd(2); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(b)");
  if (ratio)
    plotKERatio("Cu",ene,angle,first,logy,ymin,ymax,particle,beam,error,leg2,dir,dird,markf);
  else
    plotKE("Cu",ene,angle,first,logy,ymin,ymax,particle,beam,leg2,dir,dird,markf);
  myc->cd(3); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(c)");
  if (ratio)
    plotKERatio("Pb",ene,angle,first,logy,ymin,ymax,particle,beam,error,leg, dir,dird,markf);
  else
    plotKE("Pb",ene,angle,first,logy,ymin,ymax,particle,beam,leg, dir,dird,markf);
  myc->cd(4); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(d)");
  if (ratio)
    plotKERatio("U", ene,angle,first,logy,ymin,ymax,particle,beam,error,leg, dir,dird,markf);
  else
    plotKE("U", ene,angle,first,logy,ymin,ymax,particle,beam,leg, dir,dird,markf);

  char anglx[6], fname[160];
  int nx = 0;
  for (int i=0; i<6; i++) {
    if (angle[i] != ' ') { anglx[nx] = angle[i]; nx++;}
  }
  if (save != 0) {
    std::string tag=".gif";
    if (ratio) {
      if (save > 0) tag = "R.eps";
      else          tag = "R.gif";
    } else {
      if (save > 0) tag = ".eps";
    }
    sprintf (fname, "%sCCuPbUto%sat%sGeV%sdeg%s", beam, particle, ene, anglx, tag.c_str());
    myc->SaveAs(fname);
  }
}

void plotKEn(char ene[6], int first=0, int logy=0, int save=0, double ymin=-1.,
	     double ymax=-1., char beam[8]="proton", bool ratio=false, 
	     bool error=true, int leg1=1, int leg2=1,char dir[20]=".", 
	     char dird[40]=".", char mark=' ') {

  setStyle();  
  TCanvas *myc = new TCanvas("myc","",800,600); myc->Divide(2,2);

  int leg=leg1; if (leg2 < leg) leg=leg2;
  char markf[4]=" ";
  myc->cd(1); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(a)");
  if (ratio)
    plotKERatio("C","1.40","119.0",first,logy,ymin,ymax,"neutron",beam,error,leg1,dir,dird,markf);
  else
    plotKE("C","1.40","119.0",first,logy,ymin,ymax,"neutron",beam,leg1,dir,dird,markf);
  myc->cd(2); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(b)");
  if (ratio)
    plotKERatio("C",ene,   "119.0",first,logy,ymin,ymax,"neutron",beam,error,leg2,dir,dird,markf);
  else
    plotKE("C",ene,   "119.0",first,logy,ymin,ymax,"neutron",beam,leg2,dir,dird,markf);
  myc->cd(3); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(c)");
  if (ratio)
    plotKERatio("U","1.40","119.0",first,logy,ymin,ymax,"neutron",beam,error,leg,dir,dird,markf);
  else
    plotKE("U","1.40","119.0",first,logy,ymin,ymax,"neutron",beam,leg,dir,dird,markf);
  myc->cd(4); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(d)");
  if (ratio)
    plotKERatio("U",ene,   "119.0",first,logy,ymin,ymax,"neutron",beam,error,leg,dir,dird,markf);
  else
    plotKE("U",ene,   "119.0",first,logy,ymin,ymax,"neutron",beam,leg,dir,dird,markf);
  
  char fname[160];
  if (save != 0) {
    std::string tag=".gif";
    if (ratio) {
      if (save > 0) tag = "R.eps";
      else          tag = "R.gif";
    } else {
      if (save > 0) tag = ".eps";
    }
    sprintf (fname, "%sCUtoneutron_1%s", beam, tag.c_str());
    myc->SaveAs(fname);
  }

}

void plotKEp(char element[2], char ene[6], int first=0, int logy=0, int save=0,
	     double ymin=-1., double ymax=-1., char beam[8]="proton", 
	     bool ratio=false, bool error=true, int leg1=1, int leg2=1, 
	     char dir[20]=".", char dird[40]=".", char mark=' ') {

  setStyle();  
  TCanvas *myc = new TCanvas("myc","",800,600); myc->Divide(2,2);

  int leg=leg1; if (leg2 < leg) leg=leg2;
  char markf[4]=" ";
  myc->cd(1); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(a)");
  if (ratio)
    plotKERatio(element,"1.40"," 59.1",first,logy,ymin,ymax,"proton",beam,error,leg1,dir,dird,markf);
  else
    plotKE(element,"1.40"," 59.1",first,logy,ymin,ymax,"proton",beam,leg1,dir,dird,markf);
  myc->cd(2); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(b)");
  if (ratio)
    plotKERatio(element,ene,   " 59.1",first,logy,ymin,ymax,"proton",beam,error,leg2,dir,dird,markf);
  else
    plotKE(element,ene,   " 59.1",first,logy,ymin,ymax,"proton",beam,leg2,dir,dird,markf);
  myc->cd(3); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(c)");
  if (ratio)
    plotKERatio(element,"1.40","119.0",first,logy,ymin,ymax,"proton",beam,error,leg,dir,dird,markf);
  else
    plotKE(element,"1.40","119.0",first,logy,ymin,ymax,"proton",beam,leg,dir,dird,markf);
  myc->cd(4); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(d)");
  if (ratio)
    plotKERatio(element,ene,   "119.0",first,logy,ymin,ymax,"proton",beam,error,leg,dir,dird,markf);
  else
    plotKE(element,ene,   "119.0",first,logy,ymin,ymax,"proton",beam,leg,dir,dird,markf);

  char fname[160];
  if (save != 0) {
    std::string tag=".gif";
    if (ratio) {
      if (save > 0) tag = "R.eps";
      else          tag = "R.gif";
    } else {
      if (save > 0) tag = ".eps";
    }
    sprintf (fname, "%s%stoproton_1%s", beam, element, tag.c_str());
    myc->SaveAs(fname);
  }

}

void plotKE4(char element[2], char ene[6], int first=0, int logy=0, int save=0,
	     double ymin=-1., double ymax=-1., char particle[8]="proton", 
	     char beam[8]="proton", bool ratio=false, bool error=true, 
	     int leg1=1, int leg2=1, char dir[20]=".", char dird[40]=".", 
	     char mark=' ') {

  setStyle();  
  TCanvas *myc = new TCanvas("myc","",800,600); myc->Divide(2,2);

  int leg=leg1; if (leg2 < leg) leg=leg2;
  char markf[4]=" ";
  myc->cd(1); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(a)");
  if (ratio)
    plotKERatio(element,ene," 59.1",first,logy,ymin,ymax,particle,beam,error,leg1,dir,dird,markf);
  else
    plotKE(element,ene," 59.1",first,logy,ymin,ymax,particle,beam,leg1,dir,dird,markf);
  myc->cd(2); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(b)");
  if (ratio)
    plotKERatio(element,ene," 89.0",first,logy,ymin,ymax,particle,beam,error,leg2,dir,dird,markf);
  else
    plotKE(element,ene," 89.0",first,logy,ymin,ymax,particle,beam,leg2,dir,dird,markf);
  myc->cd(3); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(c)");
  if (ratio)
    plotKERatio(element,ene,"119.0",first,logy,ymin,ymax,particle,beam,error,leg,dir,dird,markf);
  else
    plotKE(element,ene,"119.0",first,logy,ymin,ymax,particle,beam,leg,dir,dird,markf);
  myc->cd(4); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(d)");
  if (ratio)
    plotKERatio(element,ene,"159.6",first,logy,ymin,ymax,particle,beam,error,leg,dir,dird,markf);
  else
    plotKE(element,ene,"159.6",first,logy,ymin,ymax,particle,beam,leg,dir,dird,markf);

  char fname[160];
  if (save != 0) {
    std::string tag=".gif";
    if (ratio) {
      if (save > 0) tag = "R.eps";
      else          tag = "R.gif";
    } else {
      if (save > 0) tag = ".eps";
    }
    sprintf (fname, "%s%sto%sat%sGeV_1%s", beam, element, particle, ene, tag.c_str());
    myc->SaveAs(fname);
  }

}

void plotKE1(char element[2], char ene[6], char angle[6], int first=0, 
	     int logy=0, int save=0, double ymin=-1, double ymax=-1., 
	     char particle[8]="proton", char beam[8]="proton", 
	     bool ratio=false, bool error=true, int legend=1, char dir[20]=".",
	     char dird[40]=".", char markf[4]=" ") {

  setStyle();
  TCanvas *myc = new TCanvas("myc","",500,600); myc->SetLeftMargin(0.15);
  if (logy != 0) gPad->SetLogy(1);
  if (ratio) 
    plotKERatio(element,ene,angle,first,logy,ymin,ymax,particle,beam,error,legend,dir,dird,markf);
  else
    plotKE(element,ene,angle,first,logy,ymin,ymax,particle,beam,legend,dir,dird,markf);

  char anglx[6], fname[160];
  int nx = 0;
  for (int i=0; i<6; i++) {
    if (angle[i] != ' ') { anglx[nx] = angle[i]; nx++;}
  }
  if (save != 0) {
    std::string tag=".gif";
    if (ratio) {
      if (save > 0) tag = "R.eps";
      else          tag = "R.gif";
    } else {
      if (save > 0) tag = ".eps";
    }
    sprintf (fname, "%s%sto%sat%sGeV%sdeg%s", beam, element, particle, ene, anglx, tag.c_str());
    myc->SaveAs(fname);
  }
}

void plotKE(char element[2], char ene[6], char angle[6], int firstMode=0, 
	    int logy=0, double ymin=-1, double ymax=-1., 
	    char particle[8]="proton", char beam[8]="proton", int legend=1, 
	    char dir[20]=".", char dird[40]=".", char markf[4]=" ") {

  int first = 0, models = modelsITEP;
  bool mode = true;
  if (firstMode < 0) {
    mode  = false;
    models= modelsD;
  } else {
    first = firstMode;
  }
  if (debug) std::cout << "First " << first << " Models " << models << " Mode " << mode << std::endl;

  char fname[160], list[80], hname[80], titlx[100];
  TH1F *hi[9];
  int i=0, icol=1, isty=1;
  sprintf (titlx, "Kinetic Energy of %s (GeV)", particle);
  double  ymx0=1, ymi0=100., xlow=0.06, xhigh=0.26;
  if (particle == "neutron") {xlow= 0.0; xhigh=0.20;}
  for (i=0; i<models; i++) {
    icol = colModel[i]; isty = stylModel[i];
    if (mode) {
      sprintf (list, "%s", ModelsITEP[i].c_str());  
      sprintf (fname,"%s/%s%s%s%sGeV.root", dir, beam, element, list, ene);
      sprintf (list, "%s", ModelsITEPh[i].c_str());  
    } else {
      sprintf (list, "%s", ModelNameD.c_str()); 
      sprintf (fname,"%s/%s%s%s%sGeV.root", ModelDirectory[i].c_str(), beam, element, list, ene);
    }
    sprintf (hname, "KE%s0%s%s%s%sGeV%s", particle, beam, element, list, ene, angle);

    TFile *file = new TFile(fname);
    hi[i] = (TH1F*) file->Get(hname);
    
    if (debug) std::cout << "Get " << hname << " from " << fname <<" as " << hi[i] <<"\n";
            
    if (hi[i] != 0) {
      int nx = hi[i]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double xx = hi[i]->GetBinCenter(k);
	double yy = hi[i]->GetBinContent(k);
	if (xx > xlow && xx < xhigh) {
	  if (yy > ymx0) ymx0 = yy;
	  if (yy < ymi0 && yy > 0) ymi0 = yy;
	}
      }
      hi[i]->GetXaxis()->SetRangeUser(xlow, xhigh); hi[i]->SetTitle("");
      hi[i]->GetXaxis()->SetTitle(titlx);
      hi[i]->SetLineStyle(isty);  hi[i]->SetLineWidth(2); hi[i]->SetLineColor(icol);
      hi[i]->GetXaxis()->SetLabelSize(0.035);hi[i]->GetYaxis()->SetLabelSize(0.035);
    }
    //    file->Close();
  }
  
  char anglx[6];
  int nx = 0;
  for (i=0; i<6; i++) {
    if (angle[i] != ' ') { anglx[nx] = angle[i]; nx++;}
  }
  sprintf (fname, "%s/itep/%s/%s/%s%sGeV%sdeg.dat", dird, beam, particle, element, ene, anglx);
  if (debug) std::cout << "Reads data from file " << fname << "\n";
  ifstream infile;
  infile.open(fname);
  
  int     q1;
  float   m1, r1, x1[30], y1[30], stater1[30], syser1[30];
  infile >> m1 >> r1 >> q1;
  for (i=0; i<q1; i++) {
    infile >> x1[i] >> y1[i] >> stater1[i] >> syser1[i];
    syser1[i] *= y1[i];
    double err = sqrt(syser1[i]*syser1[i]+stater1[i]*stater1[i]);
    stater1[i] = err;
    if (y1[i]+stater1[i] > ymx0) ymx0 = y1[i]+stater1[i];    
    if (y1[i]-stater1[i] < ymi0 && y1[i]-stater1[i] > 0) ymi0=y1[i]-stater1[i];
    if (debug) std::cout << i << " " << x1[i] << " " << y1[i] << " " << stater1[i] << "\n";
  }
  TGraph*  gr1=0;
  if (q1 > 0) {
    gr1 = new TGraphErrors(q1,x1,y1,0,stater1);
    gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
    gr1->SetMarkerSize(1.6);
  }

  if (logy == 0) {ymx0 *= 1.5; ymi0 *= 0.8;}
  else           {ymx0 *=10.0; ymi0 *= 0.2; }
  if (ymin > 0) ymi0 = ymin;
  if (ymax > 0) ymx0 = ymax;
  for (i = 0; i<models; i++) {
    if (debug) std::cout << "Model " << i << " " << hi[i] << " " << ymi0 << " " << ymx0 << "\n";
    if (hi[i] != 0) hi[i]->GetYaxis()->SetRangeUser(ymi0,ymx0);
  }

  hi[first]->GetYaxis()->SetTitleOffset(1.6);
  hi[first]->Draw();
  for (i=0; i<models; i++) {
    if (i != first && hi[i] != 0) hi[i]->Draw("same");
  }
  if (gr1) gr1->Draw("p");

  TLegend *leg1;
  if (legend < 0) {
    leg1 = new TLegend(0.60,0.55,0.90,0.90);
  } else {
    if (markf == " ") leg1 = new TLegend(0.42,0.55,0.90,0.90);
    else              leg1 = new TLegend(0.38,0.70,0.90,0.90);
  }
  for (i=0; i<models; i++) {
    if (hi[i] != 0) {
      if (mode) sprintf (list, "%s", ModelNamesI[i].c_str()); 
      else      sprintf (list, "%s", ModelNamesD[i].c_str()); 
      leg1->AddEntry(hi[i],list,"F");
    }
  }
  char header[120], beamx[8], partx[2];
  if      (beam == "piplus")  sprintf (beamx, "#pi^{+}");
  else if (beam == "piminus") sprintf (beamx, "#pi^{-}");
  else                        sprintf (beamx, "p");
  if      (particle == "neutron") sprintf (partx, "n");
  else                            sprintf (partx, "p");
  if (legend < 0) {
    sprintf (header,"%s+A #rightarrow %s+X", beamx, partx);
  } else {
    if (markf == " ") 
      sprintf (header,"%s+%s #rightarrow %s+X at %s GeV (#theta = %s^{o})", beamx, element, partx, ene, angle);
    else 
      sprintf (header,"%s %s+%s #rightarrow %s+X at %s GeV (#theta = %s^{o})", markf, beamx, element, partx, ene, angle);
  }
  leg1->SetHeader(header); leg1->SetFillColor(0);
  leg1->SetTextSize(0.04);
  if (legend != 0) leg1->Draw("same");
  if (debug) std::cout << "End\n";
}

void plotKERatio(char element[2], char ene[6], char angle[6], int firstMode=0, 
		 int logy=0, double ymin=-1, double ymax=-1., 
		 char particle[8]="proton", char beam[8]="proton", 
		 bool error=true, int legend=1, char dir[20]=".",
		 char dird[40]=".", char markf[4]=" ") {

  //Decide the mode
  int first = 0, models = modelsITEP;
  bool mode = true;
  if (firstMode < 0) {
    mode  = false;
    models= modelsD;
  } else {
    first = firstMode;
  }
  if (debug) std::cout << "First " << first << " Models " << models << " Mode " << mode << std::endl;

  // First open the data file
  char anglx[6], fname[160];
  int nx = 0;
  for (int i=0; i<6; i++) {
    if (angle[i] != ' ') { anglx[nx] = angle[i]; nx++;}
  }
  sprintf (fname, "%s/itep/%s/%s/%s%sGeV%sdeg.dat", dird, beam, particle, element, ene, anglx);
  if (debug) std::cout << "Reads data from file " << fname << "\n";
  ifstream infile;
  infile.open(fname);

  // Read contents of the data file
  int     q1;
  float   m1, r1, x1[30], y1[30], er1[30], y2[30], er2[30], staterr, syserr;
  infile >> m1 >> r1 >> q1;
  for (i=0; i<q1; i++) {
    infile >> x1[i] >> y1[i] >> staterr >> syserr;
    syserr *= y1[i];
    er1[i]  = sqrt(syserr*syserr+staterr*staterr);
    y2[i]   = 1.;
    er2[i]  = er1[i]/y1[i];
    if (debug) std::cout << i << " " << x1[i] << " " << y1[i] << " " << er1[i] << " " << er2[i] << "\n";
  }

  char list[80], hname[80], titlx[100];
  TGraphErrors *gr[9], *gref;
  int icol=1, ityp=20;
  sprintf (titlx, "Kinetic Energy of %s (GeV)", particle);
  double  ymx0=0.1, ymi0=100., xlow=0.06, xhigh=0.26;
  if (particle == "neutron") {xlow= 0.0; xhigh=0.20;}
  for (int i=0; i<models; i++) {
    icol = colModel[i]; ityp = symbModel[i];
    if (mode) {
      sprintf (list, "%s", ModelsITEP[i].c_str());  
      sprintf (fname,"%s/%s%s%s%sGeV.root", dir, beam, element, list, ene);
      sprintf (list, "%s", ModelsITEPh[i].c_str());  
    } else {
      sprintf (list, "%s", ModelNameD.c_str()); 
      sprintf (fname,"%s/%s%s%s%sGeV.root", ModelDirectory[i].c_str(), beam, element, list, ene);
    }
    sprintf (hname, "KE%s0%s%s%s%sGeV%s", particle, beam, element, list, ene, angle);

    TFile *file = new TFile(fname);
    TH1F *hi = (TH1F*) file->Get(hname);
    if (debug) std::cout << "Get " << hname << " from " << fname <<" as " << hi <<"\n";
            
    if (hi != 0 && q1 > 0) {
      float xx[30], dx[30], rat[30], drt[30];
      int   nx = hi->GetNbinsX();
      int   np = 0;
      if (debug) std::cout << "Start with " << nx << " bins\n";
      for (int k=1; k <= nx; k++) {
	double xx1 = hi->GetBinLowEdge(k);
	double xx2 = hi->GetBinWidth(k);
	for (int j=0; j<q1; j++) {
	  if (xx1 < x1[j] && xx1+xx2 > x1[j]) {
	    double yy = hi->GetBinContent(k);
	    xx[np]    = x1[j];
	    dx[np]    = 0;
	    rat[np]   = yy/y1[j];
	    drt[np]   = er1[j]*rat[np]/y1[j];
	    if (xx[np] > xlow && xx[np] < xhigh) {
	      if (rat[np]+drt[np] > ymx0) ymx0 = rat[np]+drt[np];
	      if (rat[np]-drt[np] < ymi0) ymi0 = rat[np]-drt[np];
	    }
	    if (debug) std::cout << np << "/" << j << "/" << k << " x " << xx[np] << " (" << xx1 << ":" << xx1+xx2 << ")" << " y " << yy << "/" << y1[j] << " = " << rat[np] << " +- " << drt[np] << "\n";
	    if (!error) drt[np] = 0;
	    np++;
	    break;
	  }
	}
      }
      gr[i] = new TGraphErrors(np, xx, rat, dx, drt);
      gr[i]->GetXaxis()->SetRangeUser(xlow, xhigh); gr[i]->SetTitle("");
      gr[i]->GetXaxis()->SetTitle(titlx);
      gr[i]->GetYaxis()->SetTitle("MC/Data");
      gr[i]->SetLineStyle(stylModel[i]); gr[i]->SetLineWidth(2); 
      gr[i]->SetLineColor(icol);         gr[i]->SetMarkerColor(icol); 
      gr[i]->SetMarkerStyle(ityp);       gr[i]->SetMarkerSize(1.0); 
      gr[i]->GetXaxis()->SetLabelSize(0.035);
      gr[i]->GetYaxis()->SetLabelSize(0.035);
    } else {
      gr[i] = 0;
    }
    file->Close();
  }
  if (debug) std::cout << "Completed reading for " << models << " files\n";

  gref = new TGraphErrors(q1, xx, y2, dx, er2);
  gref->GetXaxis()->SetRangeUser(xlow, xhigh); gref->SetTitle("");
  gref->GetXaxis()->SetTitle(titlx);
  gref->GetYaxis()->SetTitle("MC/Data");
  gref->SetLineStyle(1);    gref->SetLineWidth(1); 
  gref->SetLineColor(1);    gref->SetMarkerColor(1); 
  gref->SetMarkerStyle(20); gref->SetMarkerSize(0.1); 
  gref->GetXaxis()->SetLabelSize(0.035);
  gref->GetYaxis()->SetLabelSize(0.035);

  if (logy == 0) {ymx0 *= 1.5; ymi0 *= 0.8;}
  else           {ymx0 *=10.0; ymi0 *= 0.2; }
  if (ymin > 0)   ymi0 = ymin;
  if (ymax > 0)   ymx0 = ymax;
  for (i = 0; i<models; i++) {
    if (debug) std::cout << "Model " << i << " " << gr[i] << " " << ymi0 << " " << ymx0 << "\n";
    if (gr[i] != 0) gr[i]->GetYaxis()->SetRangeUser(ymi0,ymx0);
  }

  gr[first]->GetYaxis()->SetTitleOffset(1.6);
  gr[first]->Draw("APl");
  for (i=0; i<models; i++) {
    if (i != first && gr[i] != 0) gr[i]->Draw("Pl");
  }
  if (!error) gref->Draw("P");

  TLegend *leg1;
  if (legend < 0) {
    leg1 = new TLegend(0.60,0.55,0.90,0.90);
  } else {
    if (markf == " ") leg1 = new TLegend(0.42,0.55,0.90,0.90);
    else              leg1 = new TLegend(0.38,0.70,0.90,0.90);
  }
  for (i=0; i<models; i++) {
    if (gr[i] != 0) {
      if (mode) sprintf (list, "%s", ModelNamesI[i].c_str()); 
      else      sprintf (list, "%s", ModelNamesD[i].c_str()); 
      leg1->AddEntry(gr[i],list,"lP");
    }
  }
  char header[120], beamx[8], partx[2];
  if      (beam == "piplus")  sprintf (beamx, "#pi^{+}");
  else if (beam == "piminus") sprintf (beamx, "#pi^{-}");
  else                        sprintf (beamx, "p");
  if      (particle == "neutron") sprintf (partx, "n");
  else                            sprintf (partx, "p");
  if (legend < 0) {
    sprintf (header,"%s+A #rightarrow %s+X", beamx, partx);
  } else {
    if (markf == " ") 
      sprintf (header,"%s+%s #rightarrow %s+X at %s GeV (#theta = %s^{o})", beamx, element, partx, ene, angle);
    else 
      sprintf (header,"%s %s+%s #rightarrow %s+X at %s GeV (#theta = %s^{o})", markf, beamx, element, partx, ene, angle);
  }
  leg1->SetHeader(header); leg1->SetFillColor(0);
  leg1->SetTextSize(0.04);
  if (legend != 0) leg1->Draw("same");

  xx[0]=xlow; xx[1]=xhigh; rat[0]=rat[1]=1.0;
  TGraph *gr0 = new TGraph(2, xx, rat);
  gr0->GetXaxis()->SetRangeUser(xlow, xhigh); gr0->SetTitle("");
  gr0->SetLineStyle(1);   gr0->SetLineWidth(1.4); 
  gr0->SetLineColor(1);   gr0->SetMarkerColor(1); 
  gr0->SetMarkerStyle(20);gr0->SetMarkerSize(1.6);
  gr0->Draw("l");
}

void plotCT4(char element[2], char ene[6], int first=0, int scan=1, int logy=0,
	     int save=0, char particle[8]="proton", char beam[8]="proton", 
	     int leg1=1, int leg2=1, char dir[20]=".", char dird[40]=".") {

  setStyle();
  TCanvas *myc = new TCanvas("myc","",800,600); myc->Divide(2,2);
  
  for (int i=0; i<4; i++) {
    myc->cd(i+1); // if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
    double ke = keproton[i];
    if (particle == "neutron") ke = keneutron[i];
    if (i == 0) 
      plotCT(element, ene,ke, first, scan, logy, particle,beam,leg1,dir,dird); 
    else
      plotCT(element, ene,ke, first, scan, logy, particle,beam,leg2,dir,dird); 
  }

  char fname[160];
  if (save != 0) {
    if (save > 0) sprintf (fname, "%s%sto%sat%sGeV_2.eps", beam, element, particle, ene);
    else          sprintf (fname, "%s%sto%sat%sGeV_2.gif", beam, element, particle, ene);
    myc->SaveAs(fname);
  }
}
 
void plotCT1(char element[2], char ene[6], double ke, int first=0, int scan=1,
	     int logy=0, int save=0, char particle[8]="proton", 
	     char beam[8]="proton", int legend=1, char dir[20]=".",
	     char dird[40]=".") {

  setStyle();
  TCanvas *myc = new TCanvas("myc","",800,600); myc->SetLeftMargin(0.15);
  if (logy != 0) gPad->SetLogy(1);
  plotCT(element, ene,ke, first, scan, logy, particle,beam, legend, dir,dird);

  char fname[160];
  if (save != 0) {
    if (save > 0) sprintf (fname, "%s%sto%sat%sGeV%4.2fGeV.eps", beam, element, particle, ene, ke);
    else          sprintf (fname, "%s%sto%sat%sGeV%4.2fGeV.gif", beam, element, particle, ene, ke);
    myc->SaveAs(fname);
  }
}

void plotCT(char element[2], char ene[6], double ke, int firstMode=0, 
	    int scan=1, int logy=0, char particle[8]="proton", 
	    char beam[8]="proton", int legend=1, char dir[20]=".", 
	    char dird[40]=".") {

  int first = 0, models = modelsITEP;
  bool mode = true;
  if (firstMode < 0) {
    mode  = false;
    models= modelsD;
  } else {
    first = firstMode;
  }
  if (debug) std::cout << "First " << first << " Models " << models << " Mode " << mode << std::endl;

  static double pi  = 3.1415926;
  static double deg = pi/180.; 
  if (debug) std::cout << "Scan " << scan;
  std::vector<double> angles = angleScan(scan);
  int    nn = (int)(angles.size());
  if (debug) std::cout << " gives " << nn << " angles\n";

  char fname[160], list[80], hname[80];
  TH1F *hi[9];
  int i=0, icol=1;
  double  ymx0=1, ymi0=100., xlow=-1.0, xhigh=1.0;
  for (i=0; i<models; i++) {
    icol = colModel[i];
    if (mode) {
      sprintf (list, "%s", ModelsITEP[i].c_str());  
      sprintf (fname,"%s/%s%s%s%sGeV.root", dir, beam, element, list, ene);
      sprintf (list, "%s", ModelsITEPh[i].c_str());  
    } else {
      sprintf (list, "%s", ModelNameD.c_str()); 
      sprintf (fname,"%s/%s%s%s%sGeV.root", ModelDirectory[i].c_str(), beam, element, list, ene);
    }
    sprintf (hname, "CT%s0%s%s%s%sGeV%4.2f", particle, beam, element, list, ene, ke);

    TFile *file = new TFile(fname);
    hi[i] = (TH1F*) file->Get(hname);
    if (debug) std::cout << "Get " << hname << " from " << fname <<" as " << hi[i] <<"\n";
    if (hi[i] != 0) {
      int nx = hi[i]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double xx = hi[i]->GetBinCenter(k);
	double yy = hi[i]->GetBinContent(k);
	if (xx > xlow && xx < xhigh) {
	  if (yy > ymx0)           ymx0 = yy;
	  if (yy < ymi0 && yy > 0) ymi0 = yy;
	}
      }
      if (debug) std::cout << "Y limit " << ymi0 << " " << ymx0 << " after " << i;
      hi[i]->GetXaxis()->SetRangeUser(xlow, xhigh); hi[i]->SetTitle("");
      hi[i]->SetLineStyle(1);  hi[i]->SetLineWidth(2); hi[i]->SetLineColor(icol);
    }
    //    file->Close();
  }

  int     q1, kk0=0;
  float   m1, r1, x1[30], y1[30], stater1[30], syser1[30];

  for (int kk=0; kk<nn; kk++) {
    char angle[6], anglx[6];
    sprintf (angle, "%5.1f", angles[kk]);
    int nx = 0;
    for (i=0; i<6; i++) {
      if (angle[i] != ' ') { anglx[nx] = angle[i]; nx++;}
    }
    sprintf (fname, "%s/itep/%s/%s/%s%sGeV%sdeg.dat", dird, beam, particle, element, ene, anglx);
    ifstream infile;
    infile.open(fname);
  
    infile >> m1 >> r1 >> q1;
    for (i=0; i<q1; i++) {
      float xx1, yy1, stater, syser;
      infile >> xx1 >> yy1 >> stater >> syser;
      if (xx1 > ke-0.001 && xx1 < ke+0.001) {
	x1[kk0] = cos(deg*angles[kk]);
	y1[kk0] = yy1; stater1[kk0] = stater; syser1[kk0] = syser;
	syser *= yy1;
	double err = sqrt(syser*syser+stater*stater);
	stater1[kk0] = err; 
	if (y1[kk0]+stater1[kk0] > ymx0) ymx0 = y1[kk0]+stater1[kk0];
	if (y1[kk0]-stater1[kk0] < ymi0 && y1[kk0]-stater1[kk0] > 0) ymi0 = y1[kk0]-stater1[kk0];
	kk0++;
      }
    }
    infile.close();
    if (debug) std::cout << kk << " File " << fname << " X " << x1[kk] << " Y " << y1[kk] << " DY " << stater1[kk] << "\n";
  }

  TGraph*  gr1 = 0;
  if (kk0 > 0) {
    gr1 = new TGraphErrors(kk0,x1,y1,0,stater1);
    gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
    gr1->SetMarkerSize(1.6);
  }

  if (logy == 0) {ymx0 *= 1.5; ymi0 *= 0.8;}
  else           {ymx0 *=10.0; ymi0 *= 0.2; }
  for (i = 0; i<models; i++)
    if (hi[i] != 0) hi[i]->GetYaxis()->SetRangeUser(ymi0,ymx0);
  
  hi[first]->GetYaxis()->SetTitleOffset(1.6);
  hi[first]->Draw();
  for (i=0; i<models; i++) {
    if (i != first && hi[i] != 0)  hi[i]->Draw("same");
  }
  if (gr1) gr1->Draw("p");

  TLegend *leg1;
  if (legend == 1) leg1 = new TLegend(0.15,0.70,0.62,0.90);
  else             leg1 = new TLegend(0.15,0.55,0.62,0.90);
  for (i=0; i<models; i++) {
    if (hi[i] != 0) {
      if (mode) sprintf (list, "%s", ModelNamesI[i].c_str()); 
      else      sprintf (list, "%s", ModelNamesD[i].c_str()); 
      leg1->AddEntry(hi[i],list,"F");
    }
  }
  char header[80], beamx[8], partx[2];
  if      (beam == "piplus")  sprintf (beamx, "#pi^{+}");
  else if (beam == "piminus") sprintf (beamx, "#pi^{-}");
  else                        sprintf (beamx, "p");
  if      (particle == "neutron") sprintf (partx, "n");
  else                            sprintf (partx, "p");
  sprintf (header, "%s+%s #rightarrow %s+X at %s GeV (%4.2f GeV)", beamx, element, partx, ene, ke);
  leg1->SetHeader(header); leg1->SetFillColor(0);
  leg1->SetTextSize(0.04);
  if (legend != 0) leg1->Draw();

}

void plotBE4(char element[2], int logy=0, int scan=1, int save=0, 
	     char particle[8]="proton", char beam[8]="proton", 
	     char dir[20]=".", char dird[40]=".") {

  setStyle();
  TCanvas *myc = new TCanvas("myc","",800,600); myc->Divide(2,2);

  myc->cd(1); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  plotBE(element, " 59.1", 0.11, logy, scan, particle, beam, dir, dird);
  myc->cd(2); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  plotBE(element, " 59.1", 0.21, logy, scan, particle, beam, dir, dird);
  myc->cd(3); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  plotBE(element, "119.0", 0.11, logy, scan, particle, beam, dir, dird);
  myc->cd(4); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  plotBE(element, "119.0", 0.21, logy, scan, particle, beam, dir, dird);

  char fname[160];
  if (save != 0) {
    if (save > 0) sprintf (fname, "%s%sto%s_1.eps", beam, element, particle);
    else          sprintf (fname, "%s%sto%s_1.gif", beam, element, particle);
    myc->SaveAs(fname);
  }
}

void plotBE1(char element[2], char angle[6], double ke, int logy=0, int scan=1,
	     int save=0, char particle[8]="proton", char beam[8]="proton", 
	     char dir[20]=".", char dird[40]=".") {

  setStyle();
  TCanvas *myc = new TCanvas("myc","",800,600); myc->SetLeftMargin(0.15);
  if (logy != 0) gPad->SetLogy(1);
  plotBE(element, angle, ke, logy, scan, particle, beam, dir, dird);

  char anglx[6], fname[160];
  int i=0, nx=0;
  for (i=0; i<6; i++) {
    if (angle[i] != ' ') { anglx[nx] = angle[i]; nx++;}
  }
  if (save != 0) {
    if (save>0) sprintf (fname, "%s%sto%sat%sdeg%4.2fGeV.eps", beam, element, particle, anglx, ke);
    else        sprintf (fname, "%s%sto%sat%sdeg%4.2fGeV.gif", beam, element, particle, anglx, ke);
    myc->SaveAs(fname);
  }
}

void plotBE(char element[2], char angle[6], double ke, int logy=0, int scan=1,
	    char particle[8]="proton", char beam[8]="proton", 
	    char dir[20]=".", char dird[40]=".") {

  double ene[15];
  int    nene=0;
  double      energyScan0[8] = {5.7, 6.25, 6.5, 7.0, 7.5, 8.2, 8.5, 9.0};
  double      energyScan1[7] = {6.0, 6.5, 7.0, 7.5, 8.25, 8.5, 9.0};
  double      energyScan2[10]= {1.0, 2.0, 3.0,  5.0, 6.0, 6.5,
				7.0, 7.5, 8.25, 9.0};
  if (scan == 0) {
    nene = 8;
    for (int i=0; i<nene; i++) ene[i] = energyScan0[i];
  } else if (scan <= 1) {
    nene = 7;
    for (int i=0; i<nene; i++) ene[i] = energyScan1[i];
  } else if (scan == 2) {
    nene = 10;
    for (int i=0; i<nene; i++) ene[i] = energyScan2[i];
  } else {
    nene = 8;
    for (int i=0; i<nene; i++) ene[i] = energyScan2[i];
  }
 
  char anglx[6];
  int i=0, nx=0;
  for (i=0; i<6; i++) {
    if (angle[i] != ' ') { anglx[nx] = angle[i]; nx++;}
  }

  TGraph *gr[4];
  char fname[160], list[80], hname[80];
  int j=0, icol=1, ityp=20;
  double  ymx0=1, ymi0=10000., xmi=5.0, xmx=10.0;
  if (scan > 1) { 
    xmi = 0.5;
    xmx = 9.5;
  }
  for (i=0; i<modelsITEP; i++) {
    icol = colModel[i]; ityp = symbModel[i];
    double yt[15];
    for (j=0; j<nene; j++) {
      sprintf (list, "%s", ModelsITEP[i].c_str()); 
      sprintf (fname, "%s/%s%s%s%4.2fGeV.root", dir, beam, element, list, ene[j]);
      sprintf (list, "%s", ModelsITEPh[i].c_str()); 
      sprintf (hname, "KE%s0%s%s%s%4.2fGeV%s", particle, beam, element, list, ene[j], angle);

      TFile *file = new TFile(fname);
      TH1F *hi = (TH1F*) file->Get(hname);
      if (debug) std::cout << "Get " << hname << " from " << fname <<" as " << hi <<"\n";
      int    nk=0, nx = hi->GetNbinsX();
      double yy0=0;
      for (int k=1; k <= nx; k++) {
	double xx0 = hi->GetBinCenter(k);
	if (xx0 > ke-0.01 && xx0 < ke+0.01) {
	  yy0 += hi->GetBinContent(k);
	  nk++;
	}
      }
      if (nk > 0 )                 yy0 /= nk;
      if (yy0 > ymx0)              ymx0 = yy0;
      if (yy0 < ymi0 && yy0 > 0.1) ymi0 = yy0;
      if (debug) std::cout << hname << " # " << nk << " Y " << yy0 << " min " << ymi0 << " max " << ymx0 << "\n";
      yt[j] = yy0;
      file->Close();
    }
    gr[i] = new TGraph(nene, ene, yt); gr[i]->SetMarkerSize(1.2);
    gr[i]->SetTitle(list); gr[i]->SetLineColor(icol); 
    gr[i]->SetLineStyle(i+1); gr[i]->SetLineWidth(2);
    gr[i]->SetMarkerColor(icol);  gr[i]->SetMarkerStyle(ityp); 
    gr[i]->GetXaxis()->SetTitle("Beam Energy (GeV)");
    if (debug) {
      std::cout << "Graph " << i << " with " << nene << " points\n";
      for (j=0; j<nene; j++) std::cout << j << " x " << ene[j] << " y " << yt[j] << "\n";
    }
  }

  double ye[15], dy[15];
  for (j=0; j<nene; j++) {
    sprintf (fname, "%s/itep/%s/%s/%s%3.1fGeV%sdeg.dat", dird, beam, particle, element, ene[j], anglx);
    if (debug) std::cout << "Reads data from file " << fname << "\n";
    ifstream infile;
    infile.open(fname);
  
    int     q1;
    float   m1, r1, xx, yy, stater, syser;
    infile >> m1 >> r1 >> q1;
    for (i=0; i<q1; i++) {
      infile >> xx >> yy >> stater >> syser;
      if (xx > ke-0.01 && xx < ke+0.01) {
	ye[j] = yy;
	syser *= yy;
	double err = sqrt(syser*syser+stater*stater);
	dy[j] = err;
      }
    }
    infile.close();
    if (ye[j]+dy[j] > ymx0) ymx0 = ye[j]+dy[j];
    if (ye[j]-dy[j] < ymi0 && ye[j]-dy[j] > 0) ymi0 = ye[j]-dy[j];
  }
  if (debug) {
    std::cout << "Graph Data with " << nene << " points\n";
    for (j=0; j<nene; j++) std::cout << j << " x " << ene[j] << " y " << ye[j] << " +- " << dy[j] << "\n";
  }
  TGraph*  gr1=0;
  if (nene > 0) {
    gr1 = new TGraphErrors(nene,ene,ye,0,dy);
    gr1->SetMarkerColor(1);  gr1->SetMarkerStyle(22);
    gr1->SetMarkerSize(1.6);
  }

  if (logy == 0) {
    ymx0 *= 1.8; ymi0 *= 0.8;
  } else {
    ymx0 *= 50.0; ymi0 *= 0.2;
    if (scan > 1) ymx0 *= 4;
  }
  for (i = 0; i<modelsITEP; i++) {
    gr[i]->GetYaxis()->SetRangeUser(ymi0,ymx0);
    gr[i]->GetXaxis()->SetRangeUser(xmi,xmx);
  }
  if (gr1) {
    gr1->GetXaxis()->SetRangeUser(xmi,xmx);
    gr1->GetYaxis()->SetRangeUser(ymi0,ymx0);
    gr1->GetXaxis()->SetTitle("Energy (GeV)"); 
    gr1->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})"); 

    gr1->GetYaxis()->SetTitleOffset(1.6); gr1->SetTitle("");
    gr1->Draw("ap");
    for (i=0; i<modelsITEP; i++)
      gr[i]->Draw("lp");
  
    TLegend *leg1 = new TLegend(0.35,0.60,0.90,0.90);
    for (i=0; i<modelsITEP; i++) {
      sprintf (list, "%s", ModelNamesI[i].c_str());
      leg1->AddEntry(gr[i],list,"LP");
    }
    char header[80], beamx[8], partx[2];
    if      (beam == "piplus")  sprintf (beamx, "#pi^{+}");
    else if (beam == "piminus") sprintf (beamx, "#pi^{-}");
    else                        sprintf (beamx, "p");
    if      (particle == "neutron") sprintf (partx, "n");
    else                            sprintf (partx, "p");
    sprintf (header, "%s+%s #rightarrow %s+X at (KE = %3.1f GeV, #theta = %s^{o})", beamx, element, partx, ke, angle);
    leg1->SetHeader(header); leg1->SetFillColor(0);
    leg1->SetTextSize(0.04);
    leg1->Draw();
  }
}
 
void plotMT4(char ene[6], int first=0, int logy=0, int save=0, double ymin=-1,
	     double ymax=-1., char particle[8]="piplus", char beam[8]="proton",
	     bool ratio=false, bool error=true, int leg1=1, int leg2=1, 
	     char dir[20]=".", char dird[40]=".", char mark=' ') {

  setStyle();
  TCanvas *myc = new TCanvas("myc","",800,600); myc->Divide(2,2);

  char markf[4]=" ";
  myc->cd(1); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(a)");
  if (ratio)
    plotMTRatio("Be",ene,"1.10", first,logy,ymin,ymax,particle,beam,error,leg1,dir,dird,markf);
  else
    plotMT("Be",ene,"1.10", first,logy,ymin,ymax,particle,beam,leg1,dir,dird,markf);
  myc->cd(2); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(b)");
  if (ratio)
    plotMTRatio("Be",ene,"2.30", first,logy,ymin,ymax,particle,beam,error,leg2,dir,dird,markf);
  else
    plotMT("Be",ene,"2.30", first,logy,ymin,ymax,particle,beam,leg2,dir,dird,markf);
  myc->cd(3); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(c)");
  if (ratio)
    plotMTRatio("Au",ene,"1.10", first,logy,ymin,ymax,particle,beam,error,leg2,dir,dird,markf);
  else
    plotMT("Au",ene,"1.10", first,logy,ymin,ymax,particle,beam,leg2,dir,dird,markf);
  myc->cd(4); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(d)");
  if (ratio)
    plotMTRatio("Au",ene,"2.30", first,logy,ymin,ymax,particle,beam,error,leg2,dir,dird,markf);
  else
    plotMT("Au",ene,"2.30", first,logy,ymin,ymax,particle,beam,leg2,dir,dird,markf);

  char fname[160];
  if (save != 0) {
    std::string tag=".gif";
    if (ratio) {
      if (save > 0) tag = "R.eps";
      else          tag = "R.gif";
    } else {
      if (save > 0) tag = ".eps";
    }
    sprintf (fname, "%sBeAuto%sat%sGeV%s", beam, particle, ene, tag.c_str());
    myc->SaveAs(fname);
  }
}
 
void plotMT4(char element[2], char ene[6], int first=0, int logy=0, int save=0,
	     double ymin=-1, double ymax=-1., char particle[8]="piplus", 
	     char beam[8]="proton", bool ratio=false, bool error=true, 
	     int leg1=1, int leg2=1, char dir[20]=".", char dird[40]=".", 
	     char mark=' ') {

  setStyle();
  TCanvas *myc = new TCanvas("myc","",800,600); myc->Divide(2,2);

  char markf[4]=" ";
  myc->cd(1); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(a)");
  if (ratio)
    plotMTRatio(element,ene,"1.10",first,logy,ymin,ymax,particle,beam,error,leg1,dir,dird,markf);
  else
    plotMT(element,ene,"1.10",first,logy,ymin,ymax,particle,beam,leg1,dir,dird,markf);
  myc->cd(2); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(b)");
  if (ratio)
    plotMTRatio(element,ene,"1.50",first,logy,ymin,ymax,particle,beam,error,leg2,dir,dird,markf);
  else
    plotMT(element,ene,"1.50",first,logy,ymin,ymax,particle,beam,leg2,dir,dird,markf);
  myc->cd(3); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(c)");
  if (ratio)
    plotMTRatio(element,ene,"1.90",first,logy,ymin,ymax,particle,beam,error,leg2,dir,dird,markf);
  else
    plotMT(element,ene,"1.90",first,logy,ymin,ymax,particle,beam,leg2,dir,dird,markf);
  myc->cd(4); if (logy != 0) gPad->SetLogy(1); gPad->SetLeftMargin(0.15);
  if (mark == 'y') sprintf(markf, "(d)");
  if (ratio)
    plotMTRatio(element,ene,"2.30",first,logy,ymin,ymax,particle,beam,error,leg2,dir,dird,markf);
  else
    plotMT(element,ene,"2.30",first,logy,ymin,ymax,particle,beam,leg2,dir,dird,markf);

  char fname[160];
  if (save != 0) {
    std::string tag=".gif";
    if (ratio) {
      if (save > 0) tag = "R.eps";
      else          tag = "R.gif";
    } else {
      if (save > 0) tag = ".eps";
    }
    sprintf (fname, "%s%sto%sat%sGeV%s", beam, element, particle, ene, tag.c_str());
    myc->SaveAs(fname);
  }
}
 
void plotMT1(char element[2], char ene[6], char rapid[6], int first=0, 
	     int logy=0, int save=0, double ymin=-1, double ymax=-1., 
	     char particle[8]="piplus", char beam[8]="proton",bool ratio=false,
	     bool error=true, int legend=1, char dir[20]=".",char dird[40]=".",
	     char markf[4]=" ") {

  setStyle();
  TCanvas *myc = new TCanvas("myc","",500,600); myc->SetLeftMargin(0.15);
  if (logy != 0) gPad->SetLogy(1);
  if (ratio)
    plotMTRatio(element,ene,rapid,first,logy,ymin,ymax,particle,beam,error,legend,dir,dird,markf);
  else
    plotMT(element,ene,rapid,first,logy,ymin,ymax,particle,beam,legend,dir,dird,markf);

  char fname[160];
  if (save != 0) {
    std::string tag=".gif";
    if (ratio) {
      if (save > 0) tag = "R.eps";
      else          tag = "R.gif";
    } else {
      if (save > 0) tag = ".eps";
    }
    sprintf (fname, "%s%sto%sat%sGeVY%s%s", beam, element, particle, ene, rapid, tag.c_str());
    myc->SaveAs(fname);
  }
}

void plotMT(char element[2], char ene[6], char rapid[6], int firstMode=0, 
	    int logy=0, double ymin=-1, double ymax=-1., 
	    char particle[8]="piplus", char beam[8]="proton", int legend=0, 
	    char dir[20]=".", char dird[40]=".", char markf[4]=" ") {

  int first = 0, models = modelsBNL;
  bool mode = true;
  if (firstMode < 0) {
    mode  = false;
    models= modelsD;
  } else {
    first = firstMode;
  }
  if (debug) std::cout << "First " << first << " Models " << models << " Mode " << mode << std::endl;

  char fname[160], list[80], hname[80], titlx[100], sym[8];
  TH1F *hi[9];
  int i=0, icol=1;
  if      (particle=="piminus") sprintf(sym, "#pi^{-}");
  else if (particle=="piplus")  sprintf(sym, "#pi^{+}");
  else if (particle=="kminus")  sprintf(sym, "K^{-}");
  else if (particle=="kplus")   sprintf(sym, "K^{+}");
  else                          sprintf(sym, "p");
  sprintf (titlx, "Reduced m_{T} (GeV)");
  double  ymx0=1, ymi0=100., xlow=0.1, xhigh=1.6;
  for (i=0; i<models; i++) {
    icol = colModel[i];
    if (mode) {
      sprintf (list, "%s", ModelsBNL[i].c_str());
      sprintf (fname,"%s/%s%s%s%sGeV.root", dir, beam, element, list, ene);
      sprintf (list, "%s", ModelsBNLh[i].c_str());
    } else {
      sprintf (list, "%s", ModelNameD.c_str()); 
      sprintf (fname,"%s/%s%s%s%sGeV.root", ModelDirectory[i].c_str(), beam, element, list, ene);
    }
    sprintf (hname, "MT%s0%s%s%s%sGeV%s", particle, beam, element, list, ene, rapid);

    TFile *file = new TFile(fname);
    hi[i] = (TH1F*) file->Get(hname);
    if (debug) std::cout << "Get " << hname << " from " << fname <<" as " << hi[i] <<"\n";

    if (hi[i] != 0) {
      int nx = hi[i]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double xx = hi[i]->GetBinCenter(k);
	double yy = hi[i]->GetBinContent(k);
	if (xx > xlow && xx < xhigh) {
	  if (yy > ymx0) ymx0 = yy;
	  if (yy < ymi0 && yy > 0) ymi0 = yy;
	}
      }
      if (debug) std::cout << "ylimit " << ymi0 << ":" << ymx0 << "\n";
      hi[i]->GetXaxis()->SetRangeUser(xlow, xhigh); hi[i]->SetTitle("");
      hi[i]->GetXaxis()->SetTitle(titlx);
      hi[i]->SetLineStyle(1);  hi[i]->SetLineWidth(2); hi[i]->SetLineColor(icol);
    }
    //    file->Close();
  }

  sprintf (fname, "%s/bnl802/%s/%s/%s%sGeVRap%s.dat", dird, beam, particle, element, ene, rapid);
  if (debug) std::cout << "Reads data from file " << fname << "\n";
  ifstream infile;
  infile.open(fname);
  int     q1;
  float   ym1, ym2, sys, x1[50], y1[50], stater1[50], syser1[50];
  infile >> q1 >> ym1 >> ym2 >> sys;
  for (i=0; i<q1; i++) {
    infile >> x1[i] >> y1[i] >> stater1[i];
    syser1[i] = sys*y1[i];
    double err = sqrt(syser1[i]*syser1[i]+stater1[i]*stater1[i]);
    stater1[i] = err;
    if (y1[i]+stater1[i] > ymx0) ymx0 = y1[i]+stater1[i];    
    if (y1[i]-stater1[i] < ymi0 && y1[i]-stater1[i] > 0) ymi0=y1[i]-stater1[i];
    if (debug) std::cout << i << " " << x1[i] << " " << y1[i] << " " << stater1[i] << "\n";
  }
  TGraph*  gr1=0;
  if (q1 > 0) {
    gr1 = new TGraphErrors(q1,x1,y1,0,stater1);
    gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
    gr1->SetMarkerSize(1.6);
  }

  if (logy == 0) {ymx0 *= 1.5; ymi0 *= 0.8;}
  else           {ymx0 *=10.0; ymi0 *= 0.2; }
  if (ymin > 0) ymi0 = ymin;
  if (ymax > 0) ymx0 = ymax;
  for (i = 0; i<models; i++) {
    if (hi[i] != 0) {
      if (debug) std::cout << "Model " << i << " " << hi[i] << " " << ymi0 << " " << ymx0 << "\n";
      hi[i]->GetYaxis()->SetRangeUser(ymi0,ymx0);
    }
  }

  hi[first]->GetYaxis()->SetTitleOffset(1.1);
  hi[first]->Draw();
  for (i=0; i<models; i++) {
    if (i != first && hi[i] != 0) hi[i]->Draw("same");
  }
  if (gr1) gr1->Draw("p");

  TLegend *leg1;
  if (legend < 0) {
    leg1 = new TLegend(0.50,0.55,0.90,0.90);
  } else {
    if (markf == " " ) leg1 = new TLegend(0.42,0.70,0.90,0.90);
    else               leg1 = new TLegend(0.38,0.70,0.90,0.90);
  }
  for (i=0; i<models; i++) {
    if (hi[i] != 0) {
      if (mode) sprintf (list, "%s", ModelNamesB[i].c_str()); 
      else      sprintf (list, "%s", ModelNamesD[i].c_str()); 
      leg1->AddEntry(hi[i],list,"F");
    }
  }
  char header[120], beamx[8], partx[2];
  if      (beam == "piplus")  sprintf (beamx, "#pi^{+}");
  else if (beam == "piminus") sprintf (beamx, "#pi^{-}");
  else                        sprintf (beamx, "p");
  if (legend < 0) {
    sprintf (header,"%s+%s #rightarrow %s+X at %s GeV", beamx, element, sym, ene);
  } else {
    if (markf == " ")
      sprintf (header,"%s+%s #rightarrow %s+X at %s GeV (y = %s)", beamx, element, sym, ene, rapid);
    else
      sprintf (header,"%s %s+%s #rightarrow %s+X at %s GeV (y = %s)", markf, beamx, element, sym, ene, rapid);
  }
  leg1->SetHeader(header); leg1->SetFillColor(0);
  leg1->SetTextSize(0.04);
  if (legend != 0) leg1->Draw("same");

  if (debug) std::cout << "End\n";
}

void plotMTRatio(char element[2], char ene[6], char rapid[6], int firstMode=0, 
		 int logy=0, double ymin=-1, double ymax=-1., 
		 char particle[8]="piplus", char beam[8]="proton", 
		 bool error=true, int legend=0, char dir[20]=".", 
		 char dird[40]=".", char markf[4]=" ") {

  int first = 0, models = modelsBNL;
  bool mode = true;
  if (firstMode < 0) {
    mode  = false;
    models= modelsD;
  } else {
    first = firstMode;
  }
  if (debug) std::cout << "First " << first << " Models " << models << " Mode " << mode << std::endl;

  char titlx[100], sym[8];
  int  i=0, icol=1, ityp=20;
  if      (particle=="piminus") sprintf(sym, "#pi^{-}");
  else if (particle=="piplus")  sprintf(sym, "#pi^{+}");
  else if (particle=="kminus")  sprintf(sym, "K^{-}");
  else if (particle=="kplus")   sprintf(sym, "K^{+}");
  else                          sprintf(sym, "p");
  sprintf (titlx, "Reduced m_{T} (GeV)");

  //Read in the data files
  char fname[160];
  sprintf (fname, "%s/bnl802/%s/%s/%s%sGeVRap%s.dat", dird, beam, particle, element, ene, rapid);
  if (debug) std::cout << "Reads data from file " << fname << "\n";
  ifstream infile;
  infile.open(fname);
  int     q1;
  float   ym1, ym2, sys, x1[50], y1[50], y2[50], er1[50], er2[50], staterr, syserr;
  infile >> q1 >> ym1 >> ym2 >> sys;
  for (i=0; i<q1; i++) {
    infile >> x1[i] >> y1[i] >> staterr;
    syserr = sys*y1[i];
    er1[i] = sqrt(syserr*syserr+staterr*staterr);
    y2[i]  = 1.;
    er2[i] = er1[i]/y1[i];
    if (debug) std::cout << i << " " << x1[i] << " " << y1[i] << " " << er1[i] << " " << er2[i] << "\n";
  }

  char          list[80], hname[80];
  TGraphErrors *gr[9], *gref;
  double        ymx0=1, ymi0=100., xlow=0.1, xhigh=1.6;
  for (i=0; i<models; i++) {
    icol = colModel[i]; ityp = symbModel[i];
    if (mode) {
      sprintf (list, "%s", ModelsBNL[i].c_str());
      sprintf (fname, "%s/%s%s%s%sGeV.root", dir, beam, element, list, ene);
      sprintf (list, "%s", ModelsBNLh[i].c_str());
    } else {
      sprintf (list, "%s", ModelNameD.c_str()); 
      sprintf (fname,"%s/%s%s%s%sGeV.root", ModelDirectory[i].c_str(), beam, element, list, ene);
    }
    sprintf (hname, "MT%s0%s%s%s%sGeV%s", particle, beam, element, list, ene, rapid);

    TFile *file = new TFile(fname);
    TH1F  *hi   = (TH1F*) file->Get(hname);
    if (debug) std::cout << "Get " << hname << " from " << fname <<" as " << hi <<"\n";

    if (hi != 0 && q1 > 0) {
      float xx[50], dx[50], rat[50], drt[50];
      int   nx = hi->GetNbinsX();
      int   np = 0;
      if (debug) std::cout << "Start with " << nx << " bins\n";
      for (int k=1; k <= nx; k++) {
	double xx1 = hi->GetBinLowEdge(k);
	double xx2 = hi->GetBinWidth(k);
	for (int j=0; j<q1; j++) {
	  if (xx1 < x1[j] && xx1+xx2 > x1[j]) {
	    double yy = hi->GetBinContent(k);
	    xx[np]    = x1[j];
	    dx[np]    = 0;
	    rat[np]   = yy/y1[j];
	    drt[np]   = er1[j]*rat[np]/y1[j];
	    if (xx[np] > xlow && xx[np] < xhigh) {
	      if (rat[np]+drt[np] > ymx0) ymx0 = rat[np]+drt[np];
	      if (rat[np]-drt[np] < ymi0) ymi0 = rat[np]-drt[np];
	    }
	    if (debug) std::cout << np << "/" << j << "/" << k << " x " << xx[np] << " (" << xx1 << ":" << xx1+xx2 << ")" << " y " << yy << "/" << y1[j] << " = " << rat[np] << " +- " << drt[np] << "\n";
	    if (!error) drt[np] = 0;
	    np++;
	    break;
	  }
	}
      }
      gr[i] = new TGraphErrors(np, xx, rat, dx, drt);
      gr[i]->GetXaxis()->SetRangeUser(xlow, xhigh); gr[i]->SetTitle("");
      gr[i]->GetXaxis()->SetTitle(titlx);
      gr[i]->GetYaxis()->SetTitle("MC/Data");
      gr[i]->SetLineStyle(stylModel[i]); gr[i]->SetLineWidth(2); 
      gr[i]->SetLineColor(icol);         gr[i]->SetMarkerColor(icol); 
      gr[i]->SetMarkerStyle(ityp);       gr[i]->SetMarkerSize(1.0); 
    } else {
      gr[i] = 0;
    }
    file->Close();
  }
  gref = new TGraphErrors(q1, xx, y2, dx, er2);
  gref->GetXaxis()->SetRangeUser(xlow, xhigh); gref->SetTitle("");
  gref->GetXaxis()->SetTitle(titlx);
  gref->GetYaxis()->SetTitle("MC/Data");
  gref->SetLineStyle(1);    gref->SetLineWidth(1); 
  gref->SetLineColor(1);    gref->SetMarkerColor(1); 
  gref->SetMarkerStyle(20); gref->SetMarkerSize(0.1); 
  gref->GetXaxis()->SetLabelSize(0.035);
  gref->GetYaxis()->SetLabelSize(0.035);

  if (logy == 0) {ymx0 *= 1.5; ymi0 *= 0.8;}
  else           {ymx0 *=10.0; ymi0 *= 0.2; }
  if (ymin > 0) ymi0 = ymin;
  if (ymax > 0) ymx0 = ymax;
  for (i = 0; i<models; i++) {
    if (gr[i] != 0) {
      if (debug) std::cout << "Model " << i << " " << gr[i] << " " << ymi0 << " " << ymx0 << "\n";
      gr[i]->GetYaxis()->SetRangeUser(ymi0,ymx0);
    }
  }

  if (gr[first] > 0) {
    gr[first]->GetYaxis()->SetTitleOffset(1.1);
    gr[first]->Draw("APl");
    for (i=0; i<models; i++) {
      if (i != first && gr[i] != 0) gr[i]->Draw("Pl");
    }
    if (!error) gref->Draw("P");

    TLegend *leg1;
    if (legend < 0) {
      leg1 = new TLegend(0.50,0.55,0.90,0.90);
    } else {
      if (markf == " " ) leg1 = new TLegend(0.42,0.70,0.90,0.90);
      else               leg1 = new TLegend(0.38,0.70,0.90,0.90);
    }
    for (i=0; i<models; i++) {
      if (gr[i] != 0) {
	if (mode) sprintf (list, "%s", ModelNamesB[i].c_str()); 
	else      sprintf (list, "%s", ModelNamesD[i].c_str()); 
	leg1->AddEntry(gr[i],list,"lP");
      }
    }
    char header[120], beamx[8], partx[2];
    if      (beam == "piplus")  sprintf (beamx, "#pi^{+}");
    else if (beam == "piminus") sprintf (beamx, "#pi^{-}");
    else                        sprintf (beamx, "p");
    if (legend < 0) {
      sprintf (header,"%s+%s #rightarrow %s+X at %s GeV", beamx, element, sym, ene);
    } else {
      if (markf == " ")
	sprintf (header,"%s+%s #rightarrow %s+X at %s GeV (y = %s)", beamx, element, sym, ene, rapid);
      else
	sprintf (header,"%s %s+%s #rightarrow %s+X at %s GeV (y = %s)", markf, beamx, element, sym, ene, rapid);
    }
    leg1->SetHeader(header); leg1->SetFillColor(0);
    leg1->SetTextSize(0.04);
    if (legend != 0) leg1->Draw("same");

    xx[0]=xlow; xx[1]=xhigh; rat[0]=rat[1]=1.0;
    TGraph *gr0 = new TGraph(2, xx, rat);
    gr0->GetXaxis()->SetRangeUser(xlow, xhigh); gr0->SetTitle("");
    gr0->SetLineStyle(1);   gr0->SetLineWidth(1.4); 
    gr0->SetLineColor(1);   gr0->SetMarkerColor(1); 
    gr0->SetMarkerStyle(20);gr0->SetMarkerSize(1.6);
    gr0->Draw("l");
  }
}

void setStyle() {

  gStyle->SetCanvasBorderMode(0); gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);    gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);   gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);   gStyle->SetFrameLineWidth(1);
  gStyle->SetTitleOffset(1.6,"Y");  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(1);

}

std::vector<double> angleScan(int scan) {

  std::vector<double> tmp;
  if (scan <= 1) {
    tmp.push_back(59.1);
    tmp.push_back(89.0);
    tmp.push_back(119.0);
    tmp.push_back(159.6);
  } else {
    tmp.push_back(10.1);
    tmp.push_back(15.0);
    tmp.push_back(19.8);
    tmp.push_back(24.8);
    tmp.push_back(29.5);
    tmp.push_back(34.6);
    tmp.push_back(39.6);
    tmp.push_back(44.3);
    tmp.push_back(49.3);
    tmp.push_back(54.2);
    tmp.push_back(59.1);
    tmp.push_back(64.1);
    tmp.push_back(69.1);
    tmp.push_back(74.1);
    tmp.push_back(79.1);
    tmp.push_back(84.1);
    tmp.push_back(89.0);
    tmp.push_back(98.9);
    tmp.push_back(108.9);
    tmp.push_back(119.0);
    tmp.push_back(129.1);
    tmp.push_back(139.1);
    tmp.push_back(149.3);
    tmp.push_back(159.6);
    tmp.push_back(161.4);
    tmp.push_back(165.5);
    tmp.push_back(169.5);
    tmp.push_back(173.5);
    tmp.push_back(177.0);
  }
  if (debug) {
    std::cout << "Scan " << tmp.size() << " angular regions:\n";
    for (unsigned int i=0; i<tmp.size(); i++) {
      std::cout << tmp[i];
      if (i == tmp.size()-1) std::cout << " degrees\n";
      else                   std::cout << ", ";
    }
  }
  return tmp;
}
