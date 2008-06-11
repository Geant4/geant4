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


int main(int argc, char** argv)
{

  // Control on input

  if(argc < 3) {
    std::cout << "Input parameters are not specified! Exit" << std::endl;
    exit(1);
  }

  std::string fname = argv[1];
  std::string ext = ".out";
  std::string finName = fname + ext;

  std::string refer = argv[2];

  gROOT->Reset();

  const int n = 3000;
  double x, y;

  const int n_exp = 41;
  double x_exp[n], y_exp[n];

  TCanvas *c1 = new TCanvas("c1", "c1",6,6,800,600);
  gStyle->SetOptStat(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetFrameBorderMode(0);

  std::string str1 = "p (110 MeV) in Water, Geant4 ";
  std::string str_title = str1 + " " + refer;

  const char* histTitle = str_title.c_str();
      
  TH1  *h = new TH2F("h", histTitle,120,0,120,11,0,1.1);
  h->SetLineStyle(2);
  h->GetYaxis()->SetLabelFont(132);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetLabelFont(132);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitle("z, mm");
  h->GetYaxis()->SetTitle("dE/dx, arb. units");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->Draw();

  std::ifstream in;

  in.open(finName.c_str());
  if( !in.is_open()) {
    std::cout << "Input file<" << finName << "> does not exist! Exit" << std::endl;
    exit(1);
  }
   
  TH1 *h1 = new TH1F("h1","Energy deposition (MeV/mm/event) in water",3000,0,300);
  h1->SetLineStyle(1);
  h1->SetLineWidth(2);
  h1->SetLineColor(2);

  // Ignore first blank line
  char buffer[256];
  in.getline(buffer,256);

  for (Int_t i = 0; i < n; i++) {
    in >> x >> y;
    if (!in.good()) break;
    h1->SetBinContent(i+1, y);
  }

  Double_t maxJ, maxX, maxY;
  in >> maxJ >> maxX >> maxY;
  h1->Scale(1./maxY);
  h1->Draw("histosame");

  in.close();

  std::string fname2 = "H-110MeV-endep-EXP-norm-max.txt";
  in.open(fname2.c_str());

  if( !in.is_open()) { 
    std::cout << "Input file<" << fname2 << "> does not exist! Exit" << std::endl;
    exit(1);
  }
   
  // Ignore first blank line
  in.getline(buffer,256);

  for (Int_t i=0; i<n_exp; i++) {
    in >> x_exp[i] >> y_exp[i];
    if (!in.good()) break;
    x_exp[i] = 10.*x_exp[i];
  }

  TGraph *gr = new TGraph(n_exp,x_exp,y_exp);
  gr->SetMarkerStyle(22);
  gr->SetMarkerSize(1.2);
  gr->Draw("P");
 
  // draw the legend
  TLegend *leg = new TLegend(0.2,0.65,0.45,0.86);
  leg->SetTextFont(52);
  leg->SetTextSize(0.035);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillStyle(0);
  leg->SetMargin(0.4);
  leg->SetBorderSize(1);
  TLegendEntry *entry=leg->AddEntry(gr,"Data","p");
  entry->SetTextAlign(22);
  entry=leg->AddEntry("NULL","QBBC","l");
  entry->SetLineColor(2);
  entry->SetLineStyle(1);
  entry->SetLineWidth(2);
  entry->SetTextAlign(22);
  entry->SetTextColor(1);
  leg->Draw();

  c1->Modified();
  c1->cd();
  c1->Print("p_water.gif");

  in.close();
}
