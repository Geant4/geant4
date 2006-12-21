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


void pCineBarSum()
{

  //  std::ifstream expRead;

  std::ifstream simRead;


  TString confirm;

  Int_t    iLoss, iStat, iStatMax, nGamma;
  Double_t energyGG[200], energyLossG4[200],energyLossBar[200], energyLossHPW[200], Ebin, 
           delta, delta1, delta2, delta3, step, y, pos;
  Double_t intProb[200], colDist, sum, fact, GF, lambda, aaa;

  Double_t    totalGG[200], spectrumG4[200], spectrumBar[200], spectrumHPW[200];



  Double_t delExp[200], inelGG[200], distrG4[200], 
           errX[200], errY[200], deltaBin, sumPAI, sumExp, 
           meanSim, meanExp;

  Int_t numberOfExpPoints;

  Double_t Plab, Pmax, Pmin, Tkin, Tmax, Tmin, systUp, systDw, statUp, statDw;

  Double_t massPr = 0.938; // in GeV

  Double_t sigma[200],sigmaUp[200], sigmaDw[200], 
           kinEnergy[200],kinEnergyUp[200],kinEnergyDw[200];

  Int_t numberOfSimPoints,numberOfSimBarPoints;

  // exp. data for Xsec in mm from plot
  /*
  Double_t ptFNAL[26] = {0.12,0.13,0.14,0.15,0.18,0.2,0.3,0.4,0.5,0.6,0.7,
                         1.5,3.,5.,7.,8.,10.,12.,20.,30.,40.,50.,80.,130.,200.,300.};
  Double_t stFNAL[26] = {127.,110.,108.,103.,100.,95.,90.,92.,97.,102.,108.,
                      116.,114.,114.,110.,105.,107.,105.,105.,102.,102.,101.,102.,102.,103.,103.};
  Double_t ytFNAL[26] = {10.,10.,10.,10.,10.,10.,
                        10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
  Double_t xtFNAL[26] = {0.01,0.02,0.03,0.04,0.04,0.05,0.04,0.04,0.05,0.04,0.04,
                    0.05, 0.07, 1., 1.5, 1.5, 1., 1.5, 2., 3., 4., 5., 8., 13., 20., 30.};
  */

  Double_t piFNAL[27] = {0.13,0.15,0.18,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,
                         2.,2.5,3.5,4.,5.,6.,7.,8.,20.,30.,40.,50.,60.,200.,600.};
  Double_t siFNAL[27] = {71,67,66,65,65,72,75,79,82,82,84,82,
                         83.,82.,78.,77.,77.,78.,77.,77.,77.,76.,76.,77.,78.,79.,82.};
  Double_t yiFNAL[27] = {8.,8.,8.,8.,4.,3.,3.,4.,8.,8.,8.,8.,8.,
                         8.,8.,8.,4.,4.,8.,20.,8.,4.,8.,12.,8.,4.,8.};
  Double_t xiFNAL[27] = {0.01,0.02,0.03,0.04,0.04,0.05,0.04,0.04,0.05,0.04,0.04,0.04,
                         0.2,0.3,0.4,0.4,0.5,0.6,0.7,0.8,2.,3.,4.,5.,6.,20.,60.};


  // Barashenkov data for inelastic pC and after 0.2 GeV for nC

  Double_t piDUBNA[28] = {0.1,0.12,0.14,0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,
1.,1.5,2.,3.,5.,7.,10.,20.,50.,100.,500.,1000.};
  Double_t siDUBNA[28] = {222,216,210,211,205,208,210,208,210,214,216,228,240,248,254,257,
                          260,262,260,256,252,252,250,250,248,248,248,248};
  Double_t yiDUBNA[28] = {2,2,2,2,2,2,2,2,2,2,2,2,2,
                          2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
  Double_t xiDUBNA[28] = {0.003,0.003,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,
0.004,0.004,
0.01,0.01,0.02,0.02,0.02,0.03,0.07,0.08,0.2,0.3,0.4,0.5,0.6,2.,6.};

  // mm -> mb 

  for( i = 0; i < 27; i++ )
  {
    if (i < 12) siFNAL[i] = (siFNAL[i]-3.)*3.24; 
    else        siFNAL[i] *= 3.24;
 
    yiFNAL[i] *= 0.8;
    xiFNAL[i] *= 0.4;

  }

  // Read simulation data file


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
    totalGG[i]   = 0.0;
    inelGG[i]      = 0.0;  
  }

  simRead.open("pCtotine.gg");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead >> energyGG[i] >> totalGG[i] >> inelGG[i];  

    cout<<i<<"\t"<<energyGG[i] <<"\t"<<totalGG[i]<<"\t"<<inelGG[i]<<endl;
  }
  simRead.close();

  // Geisha data to be included
  
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

  // Barashenkov class

  simRead.open("pCine.bar");
  simRead>>numberOfSimBarPoints;

  for( i = 0; i < numberOfSimBarPoints; i++ )
  {
    simRead>>energyLossBar[i]>>spectrumBar[i];
    // cout<<i<<"\t"<<energyLossG4[i] <<"\t"<<spectrumG4[i]<<endl;
  }
  simRead.close();


  // HPW-Axen to be included
  
  for( i = 0; i < numberOfSimPoints; i++ )
  {
    energyLossHPW[i] = 0.0 ;
    spectrumHPW[i]    = 0.0   ;
  }
  simRead.open("pCinehpw.g4");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLossHPW[i]>>spectrumHPW[i];
    // cout<<i<<"\t"<<energyLossG4[i] <<"\t"<<spectrumG4[i]<<endl;
  }
  simRead.close();

  // Read experimental pC inelastic data from IHEP like file



  simRead.open("pC-inel.dat");


  simRead>>numberOfExpPoints;
  cout<<"numberOfExpPoints = "<<numberOfExpPoints<<endl;
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

    Tkin = std::sqrt(Plab*Plab + massPr*massPr) - massPr*massPr;
    Tmin = std::sqrt(Pmin*Pmin + massPr*massPr) - massPr*massPr;
    Tmax = std::sqrt(Pmax*Pmax + massPr*massPr) - massPr*massPr;

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
  gg->SetTitle("p-C inelastic cross-section");

  
  TGraph* g4gheisha = new TGraph(numberOfSimPoints, energyLossG4, spectrumG4);

  g4gheisha->SetMarkerStyle(20);
  g4gheisha->SetMarkerSize(0.3);
  g4gheisha->SetMarkerColor(kBlue);
  g4gheisha->SetLineColor(kBlue);
  g4gheisha->SetLineStyle(1);
  g4gheisha->SetLineWidth(1);
  g4gheisha->SetFillColor(42);
  g4gheisha->SetTitle("G-G total");

  TGraph* g4barash = new TGraph(numberOfSimBarPoints, energyLossBar, spectrumBar);

  g4barash->SetMarkerStyle(20);
  g4barash->SetMarkerSize(0.3);
  g4barash->SetMarkerColor(41);
  g4barash->SetLineColor(41);
  g4barash->SetLineStyle(1);
  g4barash->SetLineWidth(1);
  g4barash->SetFillColor(41);
  g4barash->SetTitle("G-G total");


  
  TGraphErrors* fluka = new TGraphErrors(26, piFNAL, siFNAL, xiFNAL, yiFNAL);
  fluka->SetMarkerStyle(20);
  fluka->SetMarkerColor(49);
  fluka->SetLineColor(49);
  // fluka->SetLineStyle(1);
  // fluka->SetLineWidth(0.001);
  fluka->SetFillColor(42);
  fluka->SetTitle("FNAL data");

  
  TGraphErrors* dubna = new TGraphErrors(28, piDUBNA, siDUBNA, xiDUBNA, yiDUBNA);
  dubna->SetMarkerStyle(21);
  dubna->SetMarkerColor(40);
  dubna->SetLineColor(40);
  //dubna->SetLineStyle(1);
  // dubna->SetLineWidth(0.001);
  dubna->SetFillColor(42);
  dubna->SetTitle("DUBNA data");


  
  TGraphErrors* g4hpw = new TGraph(numberOfSimPoints, energyLossHPW, spectrumHPW);

  g4hpw->SetMarkerStyle(20);
  g4hpw->SetMarkerColor(3);
  g4hpw->SetMarkerSize(0.3);
  g4hpw->SetLineColor(3);
  g4hpw->SetTitle("G4 elastic");

    
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
    


  TCanvas* myc = new TCanvas("myc", "p-C inelastic cross-section");

  myc->SetFillStyle(0);
  myc->SetFrameFillStyle(0);

  // myc->SetFillColor(50);
  myc->SetGrid();
  gPad->SetLogx();
  // gPad->SetLogy();

  gg->Draw("al");
  g4gheisha->Draw("l");
  g4barash->Draw("l");
  g4hpw->Draw("l");
  dubna->Draw("p");
  fluka->Draw("p");
  expihep->Draw("p");

  gg->GetXaxis()->SetTitle("Proton energy (GeV)");
  gg->GetYaxis()->SetTitle("Cross-section (mb)");

  gg->GetYaxis()->SetRangeUser(190,460);
  gg->GetXaxis()->SetLimits(0.01,1100.);

  gg->Draw();

  /*
  // add the graphs to the multigraph
  TMultiGraph* mg = 
  new TMultiGraph("mg", "p-C inelastic cross-section");

  
  mg->Add(gg); 
  mg->Add(g4gheisha); 
  mg->Add(fluka); 
  mg->Add(g4hpw); 
  mg->Add(expihep); 
  mg->Add(dubna); 


  mg->Draw("ap");

  mg->GetXaxis()->SetTitle("Proton energy (GeV)");
  mg->GetYaxis()->SetTitle("Cross-section (mb)");

  // mg->Draw("ap");
  // mg->Draw("ap");
  mg->Draw();
  */

  // fluka->Draw(apb);
  

  // gg->Draw("e1p");
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

  legend->AddEntry(fluka,"fluka lectures","lpe1");
  legend->AddEntry(dubna,"Dubna Bar(1989)","lpe1");

  legend->AddEntry(g4gheisha,"Geisha inelastic","l");
  legend->AddEntry(g4barash,"Barashenkov #sigma_{in}","l");
  legend->AddEntry(gg,"G-G inelastic","l");

  legend->AddEntry(expihep,"prod-exp data","lpe1");
  legend->AddEntry(g4hpw,"G4 HPW-Axen prod","l");

  legend->Draw();
  //  gg->Draw("e1p");


}

