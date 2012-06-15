//
//
// 14.02.06, V. Grichine, first version based on Normalisation.cc 

#include <fstream>
#include <iostream>
#include <stdlib.h>

#include <TObject.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStopwatch.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TSpectrum.h>
#include <TText.h>
#include <TLatex.h>
#include <TRandom.h>


void nAineXsec()
{

  Int_t i, j, k, numOfMaterials, iSan, nbOfElements, sanIndex, row;
  Double_t maxEnergyTransfer, kineticEnergy;
  Double_t tau, gamma, bg2, beta2, rateMass, Tmax;

  std::ifstream expRead;
  std::ifstream simRead;


  TString confirm;

  Int_t    iLoss, iStat, iStatMax, nGamma;
  Double_t energyLoss[200], Ebin, delta, delta1, delta2, delta3, step, y, pos;
  Double_t intProb[200], colDist, sum, fact, GF, lambda, aaa;
  Double_t    spectrum[200] ;

  Double_t delExp[200], distr[200], errX[200], errY[200], deltaBin, sumPAI, sumExp, 
           meanSim, meanExp;
  Int_t numberOfExpPoints;

  Double_t pFNAL[7] == {34.,80.,131.,180.,215.,240.,273.};
  Double_t sFNAL[7] == {1213.,1239.,1228.,1223.,1238.,1231.,1225.};
  Double_t yFNAL[7] == {30.,11.,7.,6.,9.,9.,11.};
  Double_t xFNAL[7] == {5.,6.,7.,8.,9.,10.,11.};


  // Read experimental data file

  numberOfExpPoints = 70;
  cout << " Number of experimental points = " <<numberOfExpPoints<<endl;
  cout<<endl;
  deltaBin = 1.538;
  // cout << " Exprimental bin energy bin =" <<deltaBin<<" keV"<<endl;

  expRead.open("nCuela.g4");

  for( i = 0; i < numberOfExpPoints; i++ )
  {
    expRead>>energyLoss[i]>>distr[i];
    //  cout<<i<<"\t"<<distr[i]<<endl;
  }
  expRead.close() ;

  // Read simulation data file

  Int_t numberOfSimPoints ;

  cout<<endl ;
  numberOfSimPoints = 70;
  cout << " Number of simulation points = " <<numberOfSimPoints<<endl;
  cout<<endl ;
  Ebin = 1.538;
  cout << " Energy bin in keV for simulation file = " <<Ebin<<" keV"<<endl;
  cout<<endl;
  // Ebin *= keV;

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    energyLoss[i] = 0.0 ;
    spectrum[i]    = 0.0   ;
  }
  simRead.open("nCuine.g4");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLoss[i]>>spectrum[i];  
    // cout<<i<<"\t"<<energyLoss[i] <<"\t"<<spectrum[i]<<endl;
  }
  simRead.close();

  // Sum to get total cross-section

  
  for( i = 0; i < numberOfExpPoints; i++ ) 
  {
    spectrum[i] += distr[i];
    cout<<i<<"\t"<<energyLoss[i] <<"\t"<<spectrum[i]<<endl;
  }  



  /*
  TCanvas *c1 = new TCanvas("c1","n-Cu total cross-section",20,20,700,500);

  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetGrid();

  gPad->SetLogx();
  gPad->SetLogy();
  */


  // TGraphErrors* gr1 = new TGraphErrors(numberOfExpPoints, energyLoss, distr, errX, errY);
  TGraphErrors* gr1 = new TGraph(numberOfExpPoints, energyLoss, distr);

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(kRed);
  gr1->SetMarkerSize(0.3);
  gr1->SetLineColor(kRed);
  gr1->SetTitle("Geant4 elastic");

  
  TGraph* gr2 = new TGraph(numberOfSimPoints, energyLoss, spectrum);

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.3);
  gr2->SetMarkerColor(kBlue);
  gr2->SetLineColor(kBlue);
  gr2->SetLineStyle(1);
  gr2->SetLineWidth(1);
  gr2->SetFillColor(42);
  gr2->SetTitle("Geant4 total");
  
  TGraphErrors* gr3 = new TGraphErrors(7, pFNAL, sFNAL, xFNAL, yFNAL);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerColor(kGreen);
  gr3->SetLineColor(kGreen);
  // gr3->SetLineStyle(1);
  // gr3->SetLineWidth(0.001);
  gr3->SetFillColor(42);
  gr3->SetTitle("FNAL data");

  // add the graphs to the multigraph

  TMultiGraph* mg = 
  new TMultiGraph("mg", "n-Cu total cross-section");

  
  mg->Add(gr1); 
  mg->Add(gr2); 
  mg->Add(gr3); 

  TCanvas* myc = new TCanvas("myc", "n-Cu total cross-section");

  // myc->SetFillColor(50);
  myc->SetGrid();
  gPad->SetLogx();
  // gPad->SetLogy();

  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("Neutron energy (GeV)");
  mg->GetYaxis()->SetTitle("Cross-section (mb)");
  // mg->Draw("ap");
  mg->Draw("p");
  mg->Draw();
  

  // gr1->Draw("e1p");
  /*
  TH1F *histS = 
  new TH1F("histS","n-Cu total cross-section",
    numberOfExpPoints,0,100);

  // histS->SetMarkerStyle(20);
  histS->SetMarkerSize(0.8);
  histS->SetLineColor(kBlue);
  histS->SetStats(0);
  histS->GetXaxis()->SetTitle("Neutron energy (GeV)");
  histS->GetYaxis()->SetTitle("Cross-section (mb)");

  for( i = 0; i < numberOfExpPoints; i++ ) histS->SetBinContent(i, spectrum[i]);

  // histS->Draw("same");
  
  c1->Update();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);
  c1->Modified();
  */
  // draw the legend

  TLegend* legend = new TLegend(0.6,0.65,0.88,0.85);

  legend->SetTextFont(72);
  legend->SetTextSize(0.03);

  legend->AddEntry(gr3,"FNAL data","lpe1");
  legend->AddEntry(gr2,"Geant4 total","l");
  legend->AddEntry(gr1,"Geant4 elastic","lpe1");

  legend->Draw();
  gr1->Draw("e1p");


}

