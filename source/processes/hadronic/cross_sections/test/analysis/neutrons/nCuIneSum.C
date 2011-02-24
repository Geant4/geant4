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


void nCuIneSum()
{

  //  std::ifstream expRead;

  std::ifstream simRead;


  TString confirm;

  Int_t    i, iLoss, iStat, iStatMax, nGamma;

  Double_t energyGG[200], energyLossG4[200],energyLossHPW[200], Ebin, 
           delta, delta1, delta2, delta3, step, y, pos;

  Double_t intProb[200], colDist, sum, fact, GF, lambda, aaa;

  Double_t    totalGG[200], spectrumG4[200], spectrumHPW[200];



  Double_t delExp[200], inelGG[200], distrG4[200], 
           errX[200], errY[200], deltaBin, sumPAI, sumExp, 
           meanSim, meanExp;

  Int_t numberOfExpPoints;

  Double_t Plab, Pmax, Pmin, Tkin, Tmax, Tmin, systUp, systDw, statUp, statDw;

  Double_t massNe = 0.939565; // in GeV

  Double_t sigma[200],sigmaUp[200], sigmaDw[200], 
           kinEnergy[200],kinEnergyUp[200],kinEnergyDw[200];

  Int_t numberOfSimPoints;


  // Barashenkov data for inelastic nCu and after 0.2 GeV for nC

  Double_t piDUBNA[28] = {0.1,0.12,0.14,0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,
                          1.,1.5,2.,3.,5.,7.,10.,20.,50.,100.,500.,1000.};
  Double_t siDUBNA[28] = {810,800,780,775,770,760,760,758,765,765,770,795,810,825,830,840,848,
                          870,870,868,840,825,810,803,795,795,795,795};
  Double_t yiDUBNA[28] = {2,2,2,2,2,2,2,2,2,2,2,2,2,
                          2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
  Double_t xiDUBNA[28] = {0.003,0.003,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,
    0.004,0.004,0.01,0.01,0.02,0.02,0.02,0.03,0.07,0.08,0.2,0.3,0.4,0.5,0.6,2.,6.};

  // mm -> mb 

  // Read simulation data file


  cout<<endl ;
  numberOfSimPoints = 70;
  cout << " Number of simulation points = " <<numberOfSimPoints<<endl;
  cout<<endl ;
  Ebin = 1.538;
  cout << " Energy bin in keV for simulation file = " <<Ebin<<" keV"<<endl;
  cout<<endl;
  // Ebin *= keV;

  simRead.open("nCuTotine.gg");

  simRead>>numberOfSimPoints;

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    energyGG[i] = 0.0;
    totalGG[i]   = 0.0;
    inelGG[i]      = 0.0;  
  }
  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead >> energyGG[i] >> totalGG[i] >> inelGG[i];  

    cout<<i<<"\t"<<energyGG[i] <<"\t"<<totalGG[i]<<"\t"<<inelGG[i]<<endl;
  }
  simRead.close();

  // Geisha data to be included

  simRead.open("nCuIne.g4");

  simRead>>numberOfSimPoints;
  
  for( i = 0; i < numberOfSimPoints; i++ )
  {
    energyLossG4[i] = 0.0 ;
    spectrumG4[i]    = 0.0   ;
  }
  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLossG4[i]>>spectrumG4[i];
    // cout<<i<<"\t"<<energyLossG4[i] <<"\t"<<spectrumG4[i]<<endl;
  }
  simRead.close();


  // HPW-Axen to be included

  simRead.open("nCuInehpw.g4");

  simRead>>numberOfSimPoints;
  
  for( i = 0; i < numberOfSimPoints; i++ )
  {
    energyLossHPW[i] = 0.0 ;
    spectrumHPW[i]    = 0.0   ;
  }
  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLossHPW[i]>>spectrumHPW[i];
    // cout<<i<<"\t"<<energyLossG4[i] <<"\t"<<spectrumG4[i]<<endl;
  }
  simRead.close();

  // Read experimental nC inelastic data from IHEP like file

  simRead.open("n-Cu-inelastic.dat");


  simRead>>numberOfExpPoints;

  cout<<"numberOfExpPoints = "<<numberOfExpPoints<<endl;
  // numberOfExpPoints = 28;

  for( i = 0; i < numberOfExpPoints; i++ )
  {
    kinEnergy[i]   = 0.0;
    kinEnergyUp[i] = 0.0;
    kinEnergyDw[i] = 0.0;
    sigma[i]       = 0.0;
    sigmaUp[i]     = 0.0;
    sigmaDw[i]     = 0.0;
  }

  for( i = 0; i < numberOfExpPoints; i++ )
  {
    simRead>>iLoss>>Plab>>Pmin>>Pmax>>sigma[i]>>statUp>>statDw>>systUp>>systDw;

    Tkin = std::sqrt(Plab*Plab + massNe*massNe) - massNe*massNe;
    Tmin = std::sqrt(Pmin*Pmin + massNe*massNe) - massNe*massNe;
    Tmax = std::sqrt(Pmax*Pmax + massNe*massNe) - massNe*massNe;

    kinEnergy[i]   = Tkin;
    kinEnergyDw[i] = Tkin - Tmin;
    kinEnergyUp[i] = Tmax - Tkin;

    sigmaDw[i]     = std::sqrt(statDw*statDw + systDw*systDw);
    sigmaUp[i]     = std::sqrt(statUp*statUp + systUp*systUp);

    cout<<i<<"\t"<<kinEnergy[i] <<"\t"<<sigma[i]<<endl;
  }
  simRead.close();

  /* 
  TCanvas *c1 = new TCanvas("c1","n-C total cross-section",20,20,700,500);

  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetGrid();

  gPad->SetLogx();
  gPad->SetLogy();
  */


  // TGraphErrors* gr1 = new TGraphErrors(numberOfExpPoints, energyLoss, distr, errX, errY);

  TGraphErrors* gg = new TGraph(numberOfSimPoints, energyGG, inelGG);

  gg->SetMarkerStyle(20);
  gg->SetMarkerColor(kRed);
  gg->SetMarkerSize(0.3);
  gg->SetLineColor(kRed);
  gg->SetTitle("n-Cu inelastic cross-section");

  
  TGraph* g4gheisha = new TGraph(numberOfSimPoints, energyLossG4, spectrumG4);

  g4gheisha->SetMarkerStyle(20);
  g4gheisha->SetMarkerSize(0.3);
  g4gheisha->SetMarkerColor(kBlue);
  g4gheisha->SetLineColor(kBlue);
  g4gheisha->SetLineStyle(1);
  g4gheisha->SetLineWidth(1);
  g4gheisha->SetFillColor(42);
  g4gheisha->SetTitle("G4 inelastic");

  /*  
  TGraphErrors* gr3 = new TGraphErrors(26, piFNAL, siFNAL, xiFNAL, yiFNAL);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerColor(kGreen);
  gr3->SetLineColor(kGreen);
  // gr3->SetLineStyle(1);
  // gr3->SetLineWidth(0.001);
  gr3->SetFillColor(42);
  gr3->SetTitle("FNAL data");
  */
  
  TGraphErrors* dubna = new TGraphErrors(28, piDUBNA, siDUBNA, xiDUBNA, yiDUBNA);
  dubna->SetMarkerStyle(20);
  dubna->SetMarkerColor(3);
  dubna->SetLineColor(3);
  //dubna->SetLineStyle(1);
  // dubna->SetLineWidth(0.001);
  dubna->SetFillColor(42);
  dubna->SetTitle("DUBNA data");


  
  TGraphErrors* g4hpw = new TGraph(numberOfSimPoints, energyLossHPW, spectrumHPW);

  // g4hpw->SetMarkerStyle(20);
  // g4hpw->SetMarkerColor(6);
  // g4hpw->SetMarkerSize(0.3);
  g4hpw->SetLineColor(6);
  g4hpw->SetTitle("hpw inelastic");

    
  TGraph* expihep = new TGraphAsymmErrors(numberOfExpPoints,kinEnergy,sigma,
                                     kinEnergyDw,kinEnergyUp,sigmaDw,sigmaUp);

  expihep->SetMarkerStyle(24);
  expihep->SetMarkerSize(0.9);
  expihep->SetMarkerColor(1);
  expihep->SetLineColor(1);
  // expihep->SetLineStyle(1);
  // expihep->SetLineWidth(1);
  expihep->SetFillColor(42);
  expihep->SetTitle("prod-exp data");
    
  TCanvas* myc = new TCanvas("myc", "n-Cu inelastic cross-section");

  myc->SetFillStyle(0);
  myc->SetFrameFillStyle(0);

  // myc->SetFillColor(50);
  myc->SetGrid();
  gPad->SetLogx();
  // gPad->SetLogy();

  gg->Draw("al");
  g4gheisha->Draw("l");
  g4hpw->Draw("l");
  dubna->Draw("lp");
  expihep->Draw("p");

  gg->GetXaxis()->SetTitle("Neutron energy (GeV)");
  gg->GetYaxis()->SetTitle("Cross-section (mb)");

  gg->GetYaxis()->SetRangeUser(550,1120);
  gg->GetXaxis()->SetLimits(0.09,1100.);

  gg->Draw();

  // expihep->Draw(ap);






  // add the graphs to the multigraph
  /*     
  TMultiGraph* mg = 
  new TMultiGraph("mg", "n-C inelastic cross-section");

   
  mg->Add(gg); 
  mg->Add(g4gheisha); 
  // mg->Add(gr3); 
  mg->Add(g4hpw); 
  mg->Add(expihep); 
  mg->Add(dubna); 
  
  mg->Draw("ap");

  mg->GetXaxis()->SetTitle("Neutron energy (GeV)");
  mg->GetYaxis()->SetTitle("Cross-section (mb)");

  // mg->Draw("ap");
  // mg->Draw("p");
  mg->Draw();
    
  // gg->Draw();
  */    
  TLegend* legend = new TLegend(0.6,0.65,0.88,0.85);

  legend->SetTextFont(72);
  legend->SetTextSize(0.03);

  // legend->AddEntry(gr3,"in-exp data","lpe1");

  legend->AddEntry(dubna,"Dubna Bar(1989)","lpe1");

  legend->AddEntry(g4gheisha,"Gheisha inelastic","l");
  legend->AddEntry(gg,"G-G inelastic","l");

  legend->AddEntry(expihep,"exp database","lpe1");
  legend->AddEntry(g4hpw,"G4 HPW-Axen ","l");

  legend->Draw();

  // gg->Draw("e1p");
  
  

}

