#include <fstream>
#include <iostream>
#include <string>
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h" 
#include "TH2F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"

using namespace std;

int main(int argc, char** argv)
{

  // Control on input

  if(argc < 3) {
    cout << "Input parameters are not specified! Exit" << endl;
    exit(1);
  }

  const int nidx = 3;
  string fnm[nidx] = {"p", "he4", "c12"};;
  string tp[nidx] = {"p", "^{4}He", "^{12}C"};
  string te[nidx] = {"110", "144.3", "100"};
  string teu[nidx] = {"MeV", "MeV/u", "MeV/u"};
  string fname2[nidx] = {"H-110MeV-endep-EXP-norm-max.txt",
			 "4He-144.3MeV-endep-EXP-M03-norm-max.txt",
			 "12C100MeVen-dep-EXP-norm-max.txt"}; 

  string fname = argv[1];
  int idx = 0;
  for (; idx < nidx; idx++) {if (fname == fnm[idx]) break;}

  string refer = argv[2];

  string finName[3];
  finName[0] = fname + "_opt0.out";
  finName[1] = fname + "_opt2.out";
  finName[2] = fname + "_opt3.out";

  string legend[3] = {"QBBC opt0", "QBBC opt2", "QBBC opt3"};

  gROOT->Reset();

  const int nbin = 3000;
  double x[nbin], y[nbin];

  //  const int n_exp = 41;
  int n_exp[nidx] = {41, 25, 79};
  int nn = n_exp[idx];
  double *x_exp = new double[nn];
  double *y_exp = new double[nn];

  double maxJ, maxX, maxY, norm;
  char buffer[256];

  TCanvas *c1 = new TCanvas("c1", "c1",6,6,800,600);
  gStyle->SetOptStat(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetFrameBorderMode(0);

  ifstream in;

  in.open(fname2[idx].c_str());

  if( !in.is_open()) { 
    cout << "Input file<" << fname2[idx] << "> does not exist! Exit" << endl;
    exit(1);
  }
   
  // Ignore first blank line
  in.getline(buffer,256);

  for (int i=0; i<nn; i++) {
    in >> x_exp[i] >> y_exp[i];
    if (!in.good()) break;
    x_exp[i] = 10.*x_exp[i];
  }

  double x_max = x_exp[nn-1]*1.01;
  string hist_title = tp[idx] + " " + te[idx] + " " + teu[idx] + " " + "in Water, Geant4  " + refer;

  cout << "Data file <" << fname2[idx] << " was red " << nn << " lines" << endl;

  TH1  *h = new TH2F("h", hist_title.c_str(),100,0,x_max,11,0,1.1);
  h->SetLineStyle(2);
  h->GetYaxis()->SetLabelFont(132);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetLabelFont(132);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitle("z, mm");
  h->GetYaxis()->SetTitle("dE/dx, arb. units");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->Draw();

  TGraph *gr = new TGraph(nn,x_exp,y_exp);
  gr->SetMarkerStyle(22);
  gr->SetMarkerSize(1.2);
  gr->Draw("P");

  TLegend *leg = new TLegend(0.2,0.65,0.45,0.86);
  leg->SetTextFont(52);
  leg->SetTextSize(0.035);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillStyle(0);
  leg->SetMargin(0.4);
  leg->SetBorderSize(1);
  leg->AddEntry(gr,"Data","p");

  in.close();
  c1->Update();

  TH1D* hh[3];
  string hhh[3] = {"h0", "h1", "h2"};
  TLegendEntry *entry;

  for (int j = 0; j < 2; j++) {
    in.open(finName[j].c_str());
    if( !in.is_open()) {
      cout << "Input file<" << finName[j] << "> does not exist! Exit" << endl;
      exit(1);
    }
    cout << "File with MC <" << finName[j] << "> is opened" << endl;
   
    hh[j] = new TH1D(hhh[j].c_str(),"",nbin,0,300);
    hh[j]->SetLineStyle(1);
    hh[j]->SetLineWidth(2);
    hh[j]->SetLineColor(j+2);

    // Ignore first blank line
    in.getline(buffer,256);

    for (int k = 0; k < nbin; k++) {
      in >> x[k] >> y[k];
      if (!in.good()) break;
      hh[j]->SetBinContent(k+1, y[k]);
    }

    in >> maxJ >> maxX >> maxY;
    cout << "maxJ= " << maxJ << " maxX= " << maxX << " maxY= " << maxY << endl;
    norm = maxY;
    hh[j]->Scale(1./norm);
    hh[j]->Draw("HISTO SAME");
  
    entry=leg->AddEntry(hh[j], legend[j].c_str(), "l");
    entry->SetLineColor(j+2);
    entry->SetLineStyle(1);
    entry->SetLineWidth(2);
    entry->SetTextColor(1);
    in.close();
  }
  leg->Draw();

  c1->Update();

  //  c1->Modified();
  //  c1->cd();

  delete [] x_exp;
  delete [] y_exp;

  string fout = fnm[idx] + "_water.gif";
  c1->Print(fout.c_str());
}
