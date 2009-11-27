//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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

  double zmax[nidx] = {120., 180., 40.};
 
  string fname = argv[1];
  int idx = 0;
  for (; idx < nidx; idx++) {if (fname == fnm[idx]) break;}

  string refer = argv[2];

  string finName[3];
  finName[0] = fname + "_opt0.out";
  finName[1] = fname + "_opt2.out";
  finName[2] = fname + "_opt3.out";

  string legend[3] = {"QBBC opt0", "QBBC opt2", "QBBC opt3"};

  const int nbin = 3000;
  double x[nbin], y[nbin];

  int n_exp[nidx] = {39, 25, 76};
  int nn = n_exp[idx];
  double *x_exp = new double[nn];
  double *y_exp = new double[nn];

  double maxJ, maxX, maxY, norm;
  char buffer[256];

  gROOT->SetStyle("Plain");
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
    //    if (!in.good() ||y_exp[i] == 0.0 ) break;
    x_exp[i] = 10.*x_exp[i];
  }

  //  double x_max = x_exp[nn-1]*1.01;
  string hist_title = tp[idx] + " " + te[idx] + " " + 
    teu[idx] + " " + "in Water, Geant4  " + refer;

  cout << "Data file <" << fname2[idx] << " was red " << nn << " lines" << endl;
  
  TH1F* h0 = gPad->DrawFrame(0.0,0.0,zmax[idx],1.2,hist_title.c_str());
  h0->GetXaxis()->SetTitle("z (mm)");
  h0->GetYaxis()->SetTitle("dose (relative unit)");
  h0->Draw("AXIS SAME");
  
  TGraph *gr = new TGraph(nn,x_exp,y_exp);
  gr->SetMarkerStyle(22);
  gr->SetMarkerSize(1.2);
  gr->Draw("P");
  //  gr->Draw("P SAME");

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

  for (int j = 0; j < 3; j++) {
    in.open(finName[j].c_str());
    if( !in.is_open()) {
      cout << "Input file<" << finName[j] << "> does not exist! Exit" << endl;
      exit(1);
    }
    cout << "File with MC <" << finName[j] << "> is opened" << endl;
   
    hh[j] = new TH1D(hhh[j].c_str(),"",nbin,0.,300.);
    hh[j]->SetLineStyle(1);
    hh[j]->SetLineWidth(2);
    hh[j]->SetLineColor(j+2);

    // Ignore first blank line
    in.getline(buffer,256);

    for (int k = 0; k < nbin; k++) {
      in >> x[k] >> y[k];
      if (!in.good()) {
        cout << "Stop reading results at k= " << k << endl;
	break;
      }
      hh[j]->SetBinContent(k+1, y[k]);
      hh[j]->SetBinError(k+1, y[k]/100.);
    }

    in >> maxJ >> maxX >> maxY;
    cout << "Histo filled N= " << nbin << " maxJ= " << maxJ 
	 << " maxX= " << maxX << " maxY= " << maxY << endl;
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

  delete [] x_exp;
  delete [] y_exp;

  string fout = "A_" + fnm[idx] + "_water.gif";
  c1->Print(fout.c_str());
}
