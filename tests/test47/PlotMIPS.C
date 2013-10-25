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
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphErrors.h"

void setStyle();
void plotPL(std::string element, std::string ene, int first=0, int logy=0, 
	    double ymin=-1, double ymax=-1., std::string particle="neutron",
	    int legend=1,  std::string dir=".", std::string dird=".");
void plotXF(std::string element, std::string ene, int first=0, int logy=0, 
	    double ymin=-1, double ymax=-1., std::string particle="neutron",
	    int legend=1,  std::string dir=".", std::string dird=".");

const int modelsMIPS=2;
std::string ModelsMIPS[2]  = {"qgsp",      "ftfp" };
std::string ModelNamesM[2] = {"QGS-Preco", "FTF-Preco"};
int         colModel[9]    = {3, 2, 1, 7, 4, 6, 9, 8, 12};
int         symbModel[9]   = {24, 29, 25, 27, 26, 23, 21, 20, 22};
int         stylModel[9]   = {1, 2, 3, 4, 5, 6, 1, 7, 8};

const int modelsD=5;
//std::string ModelNameD("ftfp"), ModelNameMD("FTF-Preco");
std::string ModelNameD("qgsp"), ModelNameMD("QGS_Preco");
//std::string ModelDirectory[5] = {"9.4.p01", "9.4.ref09", "9.5.ref02", "9.5.ref09", "9.5.ref10"};
//std::string ModelNamesD[5] = {"(9.4.p01)", "(9.4.ref09)", "(9.5.ref02)", "(9.5.ref09)", "(9.5.ref10)"};
std::string ModelDirectory[5] = {"9.4.ref09", "9.5.ref00", "9.5.ref02", "9.5.ref09", "9.5.ref10"};
std::string ModelNamesD[5] = {"(9.4.ref09)", "(9.5.ref00)", "(9.5.ref02)", "(9.5.ref09)", "(9.5.ref10)"};

bool debug = false;

void plot1(std::string element, std::string ene, int first=0, int logy=0, 
	   double ymin=-1, double ymax=-1., std::string particle="neutron",
	   int legend=1, bool plPlot=true, std::string dir=".", 
	   std::string dird=".") {

  setStyle();
  TCanvas *myc = new TCanvas("myc","",500,600); myc->SetLeftMargin(0.20);
  if (logy != 0) gPad->SetLogy(1);
  if (plPlot) {
    plotPL(element, ene, first,logy,ymin,ymax,particle,legend,dir,dird);
  } else {
    plotXF(element, ene, first,logy,ymin,ymax,particle,legend,dir,dird);
  }
}

void plotPL(std::string element, std::string ene, int firstMode, int logy, 
	    double ymin, double ymax, std::string particle, 
	    int legend,  std::string dir, std::string dird) {

  int models(modelsMIPS), first(0);
  bool mode = true;
  if (firstMode < 0) {
    mode  = false;
    models= modelsD;
  } else {
    first = firstMode;
  }
  if (debug) std::cout << "First " << first << " Models " << models <<std::endl;

  char fname[200], list[100], hname[100], titlx[100], titly[100];
  TH1F *hi[9];
  int i=0, icol=1, isty=1;
  sprintf (titlx, "Momentum [GeV/c]");
  sprintf (titly, "#frac{d#sigma}{dp} [mb/(GeV/c)]");
  double  ymx0=1, ymi0=100., xlow=18.0, xhigh=99.0;
  int     nbinx=27;
  if      (ene == "56.0")  {xlow= 12.0; nbinx=19;}
  else if (ene == "120.")  {xlow= 20.0; nbinx=36;}
  xhigh   = xlow + 3.0*nbinx;
  if (debug) std::cout << "Nbin " << nbinx << " in range " << xlow << ":" << xhigh << std::endl;

  for (i=0; i<models; i++) {
    icol = colModel[i]; isty = stylModel[i];
    if (mode) {
      sprintf (list, "%s", ModelsMIPS[i].c_str());  
      sprintf (fname,"%s/proton%s%s%sGeV.root", dir.c_str(), element.c_str(), list, ene.c_str());
    } else {
      sprintf (list, "%s", ModelNameD.c_str()); 
      sprintf (fname,"%s/proton%s%s%sGeV.root", ModelDirectory[i].c_str(), element.c_str(), list, ene.c_str());
    }
    sprintf (hname, "PL%sxproton%s%s%sGeV", particle.c_str(), element.c_str(), list, ene.c_str());

    TFile *file = new TFile(fname);
    TH1F* hix = (TH1F*) file->Get(hname);
    hi[i] = new TH1F(list, titlx, nbinx, xlow, xhigh);
    
    if (debug) std::cout << "Get " << hname << " from " << fname <<" as " << hix <<"\n";
            
    if (hix != 0) {
      int nx = hix->GetNbinsX();
      double cont[100];
      double wt = hix->GetBinWidth(1)/3.0;
      for (int k=1; k<=nbinx; k++) cont[k]=0;
      for (int k=1; k <= nx; k++) {
	double xx = hix->GetBinCenter(k);
	double yy = hix->GetBinContent(k);
	if (xx > xlow && xx < xhigh) {
	  for (int j=1; j<=nbinx; j++) {
	    double xxl = hi[i]->GetBinLowEdge(j);
	    if (xx > xxl && xx < xxl+3.) {
	      double yyl = wt*yy;
	      cont[j] += yyl;
	      break;
	    }
	  }
	}
      }
      for (int k=1; k<=nbinx; k++) {
	double yy = cont[k];
	hi[i]->SetBinContent(k,yy);
	if (yy > ymx0) ymx0 = yy;
	if (yy < ymi0 && yy > 0) ymi0 = yy;
	if (debug) std::cout << "Histo[" << i << "] Bin[" << k << "] = " << yy << std::endl;
      }
      hi[i]->GetXaxis()->SetRangeUser(xlow, xhigh); hi[i]->SetTitle("");
      hi[i]->GetXaxis()->SetTitle(titlx);
      hi[i]->GetYaxis()->SetTitle(titly);
      hi[i]->SetLineStyle(isty);  hi[i]->SetLineWidth(2); hi[i]->SetLineColor(icol);
      hi[i]->GetXaxis()->SetLabelSize(0.05);hi[i]->GetYaxis()->SetLabelSize(0.05);
      hi[i]->GetXaxis()->SetTitleSize(0.05);hi[i]->GetYaxis()->SetTitleSize(0.05);
    }
    //    file->Close();
  }

  if (debug) std::cout << "Yrange : " << ymi0 << "-" << ymx0 << std::endl;

  sprintf (fname, "%s/mips/%s/plab/proton%s%sGeV.dat", dird.c_str(), particle.c_str(), element.c_str(), ene.c_str());
  if (debug) std::cout << "Reads data from file " << fname << "\n";
  ifstream infile;
  infile.open(fname);
  
  int     q1;
  float   xx1, xx2, x1[30], y1[30], err[30];
  infile >> q1;
  for (i=0; i<q1; ++i) {
    infile >> xx1 >> xx2 >> y1[i] >> err[i];
    x1[i] = 0.5*(xx1+xx2);
    if (y1[i]+err[i] > ymx0) ymx0 = y1[i]+err[i];    
    if (y1[i]-err[i] < ymi0 && y1[i]-err[i] > 0) ymi0=y1[i]-err[i];
    if (debug) std::cout << i << " " << x1[i] << " " << y1[i] << " " << err[i] << "\n";
  }
  TGraphErrors*  gr1 = new TGraphErrors(q1,x1,y1,0,err);
  gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
  gr1->SetMarkerSize(1.6);

  if (logy == 0) {ymx0 *= 1.5; ymi0 *= 0.8;}
  else           {ymx0 *=10.0; ymi0 *= 0.2; }
  if (ymin > 0) ymi0 = ymin;
  if (ymax > 0) ymx0 = ymax;
  for (i = 0; i<models; i++) {
    if (debug) std::cout << "Model " << i << " " << hi[i] << " " << ymi0 << " " << ymx0 << "\n";
    if (hi[i] != 0) hi[i]->GetYaxis()->SetRangeUser(ymi0,ymx0);
  }

  hi[first]->GetYaxis()->SetTitleOffset(1.4);
  hi[first]->Draw();
  for (i=0; i<models; i++) {
    if (i != first && hi[i] != 0) hi[i]->Draw("same");
  }
  if (gr1) gr1->Draw("p");

  TLegend *leg1 = new TLegend(0.42,0.70,0.90,0.90);
  for (i=0; i<models; i++) {
    if (hi[i] != 0) {
      if (mode) sprintf (list, "%s", ModelNamesM[i].c_str()); 
      else      sprintf (list, "%s", ModelNamesD[i].c_str()); 
      leg1->AddEntry(hi[i],list,"F");
    }
  }
  char header[200], partx[20];
  if      (particle == "neutron") sprintf (partx, "n");
  else if (particle == "piplus")  sprintf (partx, "#pi^{+}");
  else if (particle == "piminus") sprintf (partx, "#pi^{-}");
  else if (particle == "kplus")   sprintf (partx, "K^{+}");
  else if (particle == "kminus")  sprintf (partx, "K^{-}");
  else                            sprintf (partx, "p");
  sprintf (header,"p+%s #rightarrow %s+X at %s GeV/c", element.c_str(), partx, ene.c_str());
  leg1->SetHeader(header); leg1->SetFillColor(0);
  leg1->SetTextSize(0.04);
  if (legend != 0) leg1->Draw("same");
  if (debug) std::cout << "End\n";
}

void plotXF(std::string element, std::string ene, int firstMode, int logy, 
	    double ymin, double ymax, std::string particle, 
	    int legend,  std::string dir, std::string dird) {

  int models(modelsMIPS), first(0);
  bool mode = true;
  if (firstMode < 0) {
    mode  = false;
    models= modelsD;
  } else {
    first = firstMode;
  }
  if (debug) std::cout << "First " << first << " Models " << models <<std::endl;

  char fname[200], list[100], hname[100], titlx[100], titly[100];
  TH1F *hi[9];
  int i=0, icol=1, isty=1;
  sprintf (titlx, "x_{F}");
  sprintf (titly, "E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");
  double  ymx0=1, ymi0=100., xlow=0.2, xhigh=1.0;
  int     nbinx=16;
  if (debug) std::cout << "Nbin " << nbinx << " in range " << xlow << ":" << xhigh << std::endl;

  for (i=0; i<models; i++) {
    icol = colModel[i]; isty = stylModel[i];
    if (mode) {
      sprintf (list, "%s", ModelsMIPS[i].c_str());  
      sprintf (fname,"%s/proton%s%s%sGeV.root", dir.c_str(), element.c_str(), list, ene.c_str());
    } else {
      sprintf (list, "%s", ModelNameD.c_str()); 
      sprintf (fname,"%s/proton%s%s%sGeV.root", ModelDirectory[i].c_str(), element.c_str(), list, ene.c_str());
    }
    sprintf (hname, "XW%sxproton%s%s%sGeV", particle.c_str(), element.c_str(), list, ene.c_str());

    TFile *file = new TFile(fname);
    TH1F* hix = (TH1F*) file->Get(hname);
    hi[i] = new TH1F(list, titlx, nbinx, xlow, xhigh);
    
    if (debug) std::cout << "Get " << hname << " from " << fname <<" as " << hix <<"\n";
            
    if (hix != 0) {
      int nx = hix->GetNbinsX();
      double cont[100];
      double wt = hix->GetBinWidth(1)/0.05;
      if (debug) std::cout << "Bin " << nx << " Weight " << wt << "\n";
      for (int k=1; k<=nbinx; k++) cont[k]=0;
      for (int k=1; k <= nx; k++) {
	double xx = hix->GetBinCenter(k);
	double yy = (hix->GetBinContent(k));
	if (xx > xlow && xx < xhigh) {
	  for (int j=1; j<=nbinx; j++) {
	    double xxl = hi[i]->GetBinLowEdge(j);
	    if (xx > xxl && xx < xxl+0.05) {
	      double yyl = wt*yy;
	      cont[j] += yyl;
	      break;
	    }
	  }
	}
      }
      for (int k=1; k<=nbinx; k++) {
	double yy = cont[k];
	hi[i]->SetBinContent(k,yy);
	if (yy > ymx0) ymx0 = yy;
	if (yy < ymi0 && yy > 0) ymi0 = yy;
	if (debug) std::cout << "Histo[" << i << "] Bin[" << k << "] = " << yy << std::endl;
      }
      hi[i]->GetXaxis()->SetRangeUser(xlow, xhigh); hi[i]->SetTitle("");
      hi[i]->GetXaxis()->SetTitle(titlx);
      hi[i]->GetYaxis()->SetTitle(titly);
      hi[i]->SetLineStyle(isty);  hi[i]->SetLineWidth(2); hi[i]->SetLineColor(icol);
      hi[i]->GetXaxis()->SetLabelSize(0.05);hi[i]->GetYaxis()->SetLabelSize(0.05);
      hi[i]->GetXaxis()->SetTitleSize(0.05);hi[i]->GetYaxis()->SetTitleSize(0.05);
    }
    //    file->Close();
  }

  if (debug) std::cout << "Yrange : " << ymi0 << "-" << ymx0 << std::endl;

  sprintf (fname, "%s/mips/%s/xf/proton%s%sGeV.dat", dird.c_str(), particle.c_str(), element.c_str(), ene.c_str());
  if (debug) std::cout << "Reads data from file " << fname << "\n";
  int     q1=0;
  float   x1[30], y1[30], err[30];

  ifstream infile;
  infile.open(fname);

  infile >> q1;
  for (i=0; i<q1; ++i) {
    float   xx1, xx2, staterr, syserr;
    infile >> xx1 >> xx2 >> y1[i] >> staterr >> syserr;
    err[i] = sqrt(staterr*staterr+syserr*syserr);
    x1[i]  = 0.5*(xx1+xx2);
    if (y1[i]+err[i] > ymx0) ymx0 = y1[i]+err[i];    
    if (y1[i]-err[i] < ymi0 && y1[i]-err[i] > 0) ymi0=y1[i]-err[i];
    if (debug) std::cout << i << " " << x1[i] << " " << y1[i] << " " << err[i] << "\n";
  }
  TGraphErrors*  gr1 = new TGraphErrors(q1,x1,y1,0,err);
  gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
  gr1->SetMarkerSize(1.6);

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

  TLegend *leg1 = new TLegend(0.42,0.70,0.90,0.90);
  for (i=0; i<models; i++) {
    if (hi[i] != 0) {
      if (mode) sprintf (list, "%s", ModelNamesM[i].c_str()); 
      else      sprintf (list, "%s", ModelNamesD[i].c_str()); 
      leg1->AddEntry(hi[i],list,"F");
    }
  }
  char header[200], partx[20];
  if      (particle == "neutron") sprintf (partx, "n");
  else if (particle == "piplus")  sprintf (partx, "#pi^{+}");
  else if (particle == "piminus") sprintf (partx, "#pi^{-}");
  else if (particle == "kplus")   sprintf (partx, "K^{+}");
  else if (particle == "kminus")  sprintf (partx, "K^{-}");
  else                            sprintf (partx, "p");
  sprintf (header,"p+%s #rightarrow %s+X at %s GeV/c", element.c_str(), partx, ene.c_str());
  leg1->SetHeader(header); leg1->SetFillColor(0);
  leg1->SetTextSize(0.04);
  if (legend != 0) leg1->Draw("same");
  if (debug) std::cout << "End\n";
}

void setStyle() {

  gStyle->SetCanvasBorderMode(0); gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);    gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);   gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);   gStyle->SetFrameLineWidth(1);
  gStyle->SetTitleOffset(2.0,"Y");  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(1.0);
  gStyle->SetLegendBorderSize(1);

}
