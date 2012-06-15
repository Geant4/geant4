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


void piCuGGbarashG4()
{

  Int_t i, j, k, numOfMaterials, iSan, nbOfElements, sanIndex, row;
  Double_t maxEnergyTransfer, kineticEnergy;
  Double_t tau, gamma, bg2, beta2, rateMass, Tmax;

  std::ifstream expRead;
  std::ifstream simRead;


  TString confirm;

  Int_t    iLoss, iStat, iStatMax, nGamma;
  Double_t energyLossGG[200],energyLossG4[200], Ebin, delta, delta1, delta2, delta3, step, y, pos;
  Double_t intProb[200], colDist, sum, fact, GF, lambda, aaa;
  Double_t totGG[200], ineGG[200], totG4[200], ineG4[200];

  Double_t delExp[200], distr[200], errX[200], errY[200], deltaBin, 
           sumPAI, sumExp, meanSim, meanExp;
  Int_t numberOfExpPoints;

  // Barashenkov calculation
  
  Double_t pDubna[10]  = {1.,2.,3.,5.,10.,20.,50.,100.,500.,1000.};

  Double_t stDubna[10] = {1365.,1250.,1185.,1128.,1070.,1035.,1010.,1010.,1010.,1010.};

  Double_t siDubna[10] = {865.,785.,735.,705.,680.,650.,630.,630.,630.,630.};

  Double_t yDubna[10]   = {30.,11.,7.,6.,9.,9.,11.,9.,9.,11.};
  Double_t xDubna[10]   = {0.1,0.2,0.3,0.4,0.5,1.,2.,5.,13.,24.};


  Double_t pFNAL[1] == {608.};
  Double_t sFNAL[1] == {1032.};
  Double_t yFNAL[1] == {180.};
  Double_t xFNAL[1] == {5.};


  // Read experimental data file

  numberOfExpPoints = 70;
  cout << " Number of experimental points = " <<numberOfExpPoints<<endl;
  cout<<endl;
  deltaBin = 1.538;
  // cout << " Exprimental bin energy bin =" <<deltaBin<<" keV"<<endl;


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
    energyLossGG[i] = 0.0 ;
    energyLossG4[i] = 0.0 ;
    totGG[i]    = 0.0   ;
    ineGG[i]    = 0.0   ;
    totG4[i]    = 0.0   ;
    ineG4[i]    = 0.0   ;
  }

  // GG simulation

  simRead.open("piCuTotine.gg");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead >> energyLossGG[i] >> totGG[i] >> ineGG[i];  
    // cout<< i << "\t" << energyLossGG[i] <<"\t"<<totGG[i] << ineGG[i] <<endl;
  }
  simRead.close();

  // Geisha simualtion

  simRead.open("piCuIne.g4");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLossG4[i]>>ineG4[i];
   
    // cout<<i<<"\t"<< energyLossG4[i] <<"\t"<< totG4[i] << ineG4 <<endl;
  }
  simRead.close();


  simRead.open("piCuEla.g4");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLossG4[i]>>totG4[i];
    totG4[i] += ineG4[i];
    // cout<<i<<"\t"<< energyLossG4[i] <<"\t"<< totG4[i] << ineG4 <<endl;
  }
  simRead.close();



  /*
  TCanvas *c1 = new TCanvas("c1","n-Cu total cross-section",20,20,700,500);

  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetGrid();

  gPad->SetLogx();
  gPad->SetLogy();
  */


  // TGraphErrors* gr1 = new TGraphErrors(numberOfExpPoints, energyLoss, distr, errX, errY);

  // gg plots

  TGraphErrors* ggtot = new TGraph(numberOfExpPoints, energyLossGG, totGG);

  ggtot->SetMarkerStyle(20);
  ggtot->SetMarkerColor(kRed);
  ggtot->SetMarkerSize(0.3);
  ggtot->SetLineColor(kRed);
  ggtot->SetTitle("GG total");

  
  TGraph* ggine = new TGraph(numberOfSimPoints, energyLossGG, ineGG);

  ggine->SetMarkerStyle(20);
  ggine->SetMarkerSize(0.3);
  ggine->SetMarkerColor(6);
  ggine->SetLineColor(6);
  ggine->SetLineStyle(1);
  ggine->SetLineWidth(1);
  ggine->SetFillColor(42);
  ggine->SetTitle("GG inelastic");
  
  // g4 plots

  TGraphErrors* g4tot = new TGraph(numberOfExpPoints, energyLossG4, totG4);

  g4tot->SetMarkerStyle(20);
  g4tot->SetMarkerColor(kBlue);
  g4tot->SetMarkerSize(0.3);
  g4tot->SetLineColor(kBlue);
  g4tot->SetTitle("Geisha total");

  
  TGraph* g4ine = new TGraph(numberOfSimPoints, energyLossG4, ineG4);

  g4ine->SetMarkerStyle(20);
  g4ine->SetMarkerSize(0.3);
  g4ine->SetMarkerColor(38);
  g4ine->SetLineColor(38);
  g4ine->SetLineStyle(1);
  g4ine->SetLineWidth(1);
  g4ine->SetFillColor(42);
  g4ine->SetTitle("Geisha inelastic");

  // Dubna calculations
  
  TGraphErrors* bartot = new TGraphErrors(10, pDubna, stDubna, xDubna, yDubna);
  bartot->SetMarkerStyle(20);
  bartot->SetMarkerColor(kGreen);
  bartot->SetLineColor(kGreen);
  // bartot->SetLineStyle(1);
  // bartot->SetLineWidth(0.001);
  bartot->SetFillColor(42);
  bartot->SetTitle("Dubna data");

  TGraphErrors* barine = new TGraphErrors(10, pDubna, siDubna, xDubna, yDubna);
  barine->SetMarkerStyle(20);
  barine->SetMarkerColor(30);
  barine->SetLineColor(30);
  // barine->SetLineStyle(1);
  // barine->SetLineWidth(0.001);
  barine->SetFillColor(42);
  barine->SetTitle("Dubna data");

  
  TGraphErrors* fnaltot = new TGraphErrors(1, pFNAL, sFNAL, xFNAL, yFNAL);
  fnaltot->SetMarkerStyle(20);
  fnaltot->SetMarkerColor(1);
  fnaltot->SetLineColor(1);
  // fnaltot->SetLineStyle(1);
  // fnaltot->SetLineWidth(0.001);
  fnaltot->SetFillColor(42);
  fnaltot->SetTitle("Dubna data");

  // add the graphs to the multigraph

  TMultiGraph* mg = 
  new TMultiGraph("mg", "#pi^{-}-Cu cross-sections");

  
  mg->Add(ggtot); 
  mg->Add(ggine); 
  mg->Add(g4tot); 
  mg->Add(g4ine); 
  mg->Add(bartot); 
  mg->Add(barine); 
  mg->Add(fnaltot); 

  TCanvas* myc = new TCanvas("myc", "#pi^{-}-Cu cross-sections");

  // myc->SetFillColor(50);
  myc->SetGrid();
  gPad->SetLogx();
  // gPad->SetLogy();

  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("#pi^{-} energy (GeV)");
  mg->GetYaxis()->SetTitle("Cross-section (mb)");
  // mg->Draw("ap");
  mg->Draw("p");
  mg->Draw();
  

  // ggtot->Draw("e1p");
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

  legend->AddEntry(fnaltot,"FNAL total","lpe1");
  legend->AddEntry(bartot,"Dubna total","lpe1");
  legend->AddEntry(barine,"Dubna inelastic","lpe1");
  legend->AddEntry(ggine,"G-G inelastic","l");
  legend->AddEntry(ggtot,"G-G total","l");
  legend->AddEntry(g4ine,"Geisha inelastic","l");
  legend->AddEntry(g4tot,"Geisha total","l");

  legend->Draw();
  gr1->Draw("e1p");


}

