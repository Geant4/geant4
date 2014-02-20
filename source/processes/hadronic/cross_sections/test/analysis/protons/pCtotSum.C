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


void pCtotSum()
{

  Int_t i, j, k, numOfMaterials, iSan, nbOfElements, sanIndex, row;
  Double_t maxEnergyTransfer, kineticEnergy;
  Double_t tau, gamma, bg2, beta2, rateMass, Tmax;

  //  std::ifstream expRead;

  std::ifstream simRead;


  TString confirm;

  Int_t    iLoss, iStat, iStatMax, nGamma;
  Double_t energyGG[200], energyLossG4[200], Ebin, delta, delta1, delta2, delta3, step, y, pos;
  Double_t intProb[200], colDist, sum, fact, GF, lambda, aaa;

  Double_t    spectrum[200], spectrumG4[200] ;

  Double_t delExp[200], distr[200], distrG4[200], errX[200], errY[200], deltaBin, sumPAI, sumExp, 
           meanSim, meanExp;

  Int_t numberOfExpPoints;

  // exp. data for Xsec in mm from plot
  // errors are approximate

  Double_t ptFNAL[26] = {0.12,0.13,0.14,0.15,0.18,0.2,0.3,0.4,0.5,0.6,0.7,
                         1.5,3.,5.,7.,8.,10.,12.,20.,30.,40.,50.,80.,130.,200.,300.};
  Double_t stFNAL[26] = {127.,110.,108.,103.,100.,95.,90.,92.,97.,102.,108.,
                      116.,114.,114.,110.,105.,107.,105.,105.,102.,102.,101.,102.,102.,103.,103.};
  Double_t ytFNAL[26] = {10.,10.,10.,10.,10.,10.,
                        10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
  Double_t xtFNAL[26] = {0.01,0.02,0.03,0.04,0.04,0.05,0.04,0.04,0.05,0.04,0.04,
                    0.05, 0.07, 1., 1.5, 1.5, 1., 1.5, 2., 3., 4., 5., 8., 13., 20., 30.};


  // Barashenkov data for total nC
  // errors are approximate

  Double_t ptDUBNA[28] = {0.1,0.12,0.14,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,
             1.,1.5,2.,3.,5.,7.,10.,20.,50.,100.,500.,1000.};
  Double_t stDUBNA[28] = {475.,396.,350.,336.,320.,303.,294.,280.,280.,286.,296.,314.,330.,344.,
           356.,360.,364.,384.,388.,384.,364.,352.,344.,330.,324.,324.,332.,332.};
  Double_t ytDUBNA[28] = {15,14,13,12,11,10,10,10,10,10,10,9,9,
                 8.,8.,8.,8.,8.,8.,8.,8.,8.,8.,8.,8.,8.,8.,8.};
  Double_t xtDUBNA[28] = {0.003,0.003,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,
             0.004,0.004,0.01,0.01,0.02,0.02,0.02,0.03,0.07,0.08,0.2,0.3,0.4,0.5,0.6,2.,6.};

  // mm -> mb 

  for( i = 0; i < 26; i++ )
  {
    stFNAL[i] *= 3.24; 
   
 
    ytFNAL[i] *= 0.8;
    xtFNAL[i] *= 0.4;

  }

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
    energyGG[i] = 0.0;
    spectrum[i]   = 0.0;
    distr[i]      = 0.0;  
  }

  simRead.open("pCtotine.gg");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead >> energyGG[i] >> spectrum[i] >> distr[i];  

    cout<<i<<"\t"<<energyGG[i] <<"\t"<<spectrum[i]<<"\t"<<distr[i]<<endl;
  }
  simRead.close();

  // Geisha data to be included
  /*

  simRead.open("pCela.g4");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLossG4[i]>>distrG4[i];
    //  cout<<i<<"\t"<<distrG4[i]<<endl;
  }
  simRead.close() ;

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    energyLossG4[i] = 0.0 ;
    spectrumG4[i]    = 0.0   ;
  }
  simRead.open("pCine.g4");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLossG4[i]>>spectrumG4[i];
    // cout<<i<<"\t"<<energyLossG4[i] <<"\t"<<spectrumG4[i]<<endl;
  }
  simRead.close();

// Sum to get total cross-section


  for( i = 0; i < numberOfSimPoints; i++ )
  {
    spectrumG4[i] += distrG4[i];
    cout<<i<<"\t"<<energyLossG4[i] <<"\t"<<spectrumG4[i]<<endl;
  }




  
  TCanvas *c1 = new TCanvas("c1","n-C total cross-section",20,20,700,500);

  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetGrid();

  gPad->SetLogx();
  gPad->SetLogy();
  */


  // TGraphErrors* gr1 = new TGraphErrors(numberOfExpPoints, energyLoss, distr, errX, errY);
  TGraphErrors* gr1 = new TGraph(numberOfSimPoints, energyGG, distr);

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(kRed);
  gr1->SetMarkerSize(0.3);
  gr1->SetLineColor(kRed);
  gr1->SetTitle("G-G elastic");

  
  TGraph* gr2 = new TGraph(numberOfSimPoints, energyGG, spectrum);

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.3);
  gr2->SetMarkerColor(kBlue);
  gr2->SetLineColor(kBlue);
  gr2->SetLineStyle(1);
  gr2->SetLineWidth(1);
  gr2->SetFillColor(42);
  gr2->SetTitle("G-G total");
  
  TGraphErrors* gr3 = new TGraphErrors(26, ptFNAL, stFNAL, xtFNAL, ytFNAL);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerColor(kGreen);
  gr3->SetLineColor(kGreen);
  // gr3->SetLineStyle(1);
  // gr3->SetLineWidth(0.001);
  gr3->SetFillColor(42);
  gr3->SetTitle("FNAL data");

  
  TGraphErrors* gr6 = new TGraphErrors(28, ptDUBNA, stDUBNA, xtDUBNA, ytDUBNA);
  gr6->SetMarkerStyle(20);
  gr6->SetMarkerColor(40);
  gr6->SetLineColor(40);
  //gr6->SetLineStyle(1);
  // gr6->SetLineWidth(0.001);
  gr6->SetFillColor(42);
  gr6->SetTitle("DUBNA data");


  /*
  TGraphErrors* gr4 = new TGraph(numberOfSimPoints, energyLossG4, distrG4);

  gr4->SetMarkerStyle(20);
  gr4->SetMarkerColor(6);
  gr4->SetMarkerSize(0.3);
  gr4->SetLineColor(6);
  gr4->SetTitle("G4 elastic");

  
  TGraph* gr5 = new TGraph(numberOfSimPoints, energyLossG4, spectrumG4);

  gr5->SetMarkerStyle(20);
  gr5->SetMarkerSize(0.3);
  gr5->SetMarkerColor(28);
  gr5->SetLineColor(28);
  gr5->SetLineStyle(1);
  gr5->SetLineWidth(1);
  gr5->SetFillColor(42);
  gr5->SetTitle("G4 total");
  */  

  // add the graphs to the multigraph

  TMultiGraph* mg = 
  new TMultiGraph("mg", "p-C total cross-section");

  
  mg->Add(gr1); 
  mg->Add(gr2); 
  mg->Add(gr3); 
  // mg->Add(gr4); 
  // mg->Add(gr5); 
  mg->Add(gr6); 

  TCanvas* myc = new TCanvas("myc", "p-C total cross-section");

  // myc->SetFillColor(50);
  myc->SetGrid();
  gPad->SetLogx();
  // gPad->SetLogy();

  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("Proton energy (GeV)");
  mg->GetYaxis()->SetTitle("Cross-section (mb)");
  // mg->Draw("ap");
  mg->Draw("p");
  mg->Draw();
  

  // gr1->Draw("e1p");
  /*
  TH1F *histS = 
  new TH1F("histS","n-C total cross-section",
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

  legend->AddEntry(gr3,"tot-FNAL data","lpe1");
  legend->AddEntry(gr6,"Dubna B-data","lpe1");

  legend->AddEntry(gr2,"G-G total","l");
  legend->AddEntry(gr1,"G-G inelastic","lpe1");

  // legend->AddEntry(gr5,"G4 total","l");
  // legend->AddEntry(gr4,"G4 elastic","lpe1");

  legend->Draw();
  gr1->Draw("e1p");


}

