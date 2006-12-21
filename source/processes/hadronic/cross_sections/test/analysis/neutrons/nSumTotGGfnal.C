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


void nSumTotGGfnal()
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
  Double_t    spectrumCu[200], spectrumBe[200], spectrumC[200], spectrumFe[200],
              spectrumW[200], spectrumPb[200], spectrumU[200] ;

  Double_t delExp[200], distr[200], errX[200], errY[200], deltaBin, sumPAI, sumExp, 
           meanSim, meanExp;

  Int_t numberOfExpPoints;

  Double_t pUFNAL[7] == {34.,80.,131.,180.,215.,240.,273.};
  Double_t sUFNAL[7] == {3402.,3410.,3399.,3362.,3353.,3365.,3297};
  Double_t yUFNAL[7] == {113.,29.,26.,32.,39.,46.,60.};
  Double_t xUFNAL[7] == {5.,6.,7.,8.,9.,10.,11.};


  Double_t pPbFNAL[7] == {34.,80.,131.,180.,215.,240.,273.};
  Double_t sPbFNAL[7] == {2973.,2986.,2981.,2951.,2959.,2926.,2919.};
  Double_t yPbFNAL[7] == {85.,25.,21.,28.,32.,32.,48};
  Double_t xPbFNAL[7] == {5.,6.,7.,8.,9.,10.,11.};


  Double_t pWFNAL[7] == {34.,80.,131.,180.,215.,240.,273.};
  Double_t sWFNAL[7] == {2840.,2804.,2786.,2751.,2746.,2748.,2720.};
  Double_t yWFNAL[7] == {46.,16.,11.,12.,16.,16.,18.};
  Double_t xWFNAL[7] == {5.,6.,7.,8.,9.,10.,11.};


  Double_t pCuFNAL[7] == {34.,80.,131.,180.,215.,240.,273.};
  Double_t sCuFNAL[7] == {1213.,1239.,1228.,1223.,1238.,1231.,1225.};
  Double_t yCuFNAL[7] == {30.,11.,7.,6.,9.,9.,11.};
  Double_t xCuFNAL[7] == {5.,6.,7.,8.,9.,10.,11.};

  Double_t pFeFNAL[7] == {34.,80.,131.,180.,215.,240.,273.};
  Double_t sFeFNAL[7] == {1100.,1122.,1110.,1110.,1112.,1113.,1107.};
  Double_t yFeFNAL[7] == {29.,11.,7.,8.,8.,8.,10.};
  Double_t xFeFNAL[7] == {5.,6.,7.,8.,9.,10.,11.};



  Double_t pCFNAL[7] == {34.,80.,131.,180.,215.,240.,273.};
  Double_t sCFNAL[7] == {331.1,331.4,329.5,331.1,333.5,331.9,328.2};
  Double_t yCFNAL[7] == {8.6,3.4,1.7,1.5,1.8,1.8,2.1};
  Double_t xCFNAL[7] == {5.,6.,7.,8.,9.,10.,11.};



  Double_t pBeFNAL[7] == {34.,80.,131.,180.,215.,240.,273.};
  Double_t sBeFNAL[7] == {263.9,269.7,266.5,271.1,273.5,270.8,273.8};
  Double_t yBeFNAL[7] == {5.7,2.8,1.3,1.2,1.3,1.3,1.6};
  Double_t xBeFNAL[7] == {5.,6.,7.,8.,9.,10.,11.};

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
    energyLoss[i] = 0.0 ;
    spectrumU[i]    = 0.0   ;
    spectrumPb[i]    = 0.0   ;
    spectrumW[i]    = 0.0   ;
    spectrumCu[i]    = 0.0   ;
    spectrumFe[i]    = 0.0   ;
    spectrumC[i]    = 0.0   ;
    spectrumBe[i]    = 0.0   ;
  }

  // U gg simulation

  simRead.open("nUtotela.gg");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead >> energyLoss[i] >> spectrumU[i] >> distr[i];

    // cout<<i<<"\t"<<energyLoss[i] <<"\t"<<spectrum[i]<<"\t"<<distr[i]<<endl;
  }
  simRead.close();


  // Pb gg simulation

  simRead.open("nPbTotela.gg");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLoss[i]>>spectrumPb[i]>>distr[i];
    // cout<<i<<"\t"<<energyLoss[i] <<"\t"<<spectrum[i]<<endl;
  }
  simRead.close();


  // W gg simulation

  simRead.open("nWtotela.gg");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead >> energyLoss[i] >> spectrumW[i] >> distr[i];

    // cout<<i<<"\t"<<energyLoss[i] <<"\t"<<spectrum[i]<<"\t"<<distr[i]<<endl;
  }
  simRead.close();

  // Cu gg simulation

  simRead.open("nCuTotela.gg");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLoss[i]>>spectrumCu[i]>>distr[i];  
    // cout<<i<<"\t"<<energyLoss[i] <<"\t"<<spectrum[i]<<endl;
  }
  simRead.close();

  // Fe gg simulation

  simRead.open("nFeTotela.gg");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLoss[i]>>spectrumFe[i]>>distr[i];
    // cout<<i<<"\t"<<energyLoss[i] <<"\t"<<spectrum[i]<<endl;
  }
  simRead.close();


  // C gg simulation

  simRead.open("nCtotela.gg");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead >> energyLoss[i] >> spectrumC[i] >> distr[i];
  }
  simRead.close();


  // Be gg simulation

  simRead.open("nBeTotela.gg");

  for( i = 0; i < numberOfSimPoints; i++ )
  {
    simRead>>energyLoss[i]>>spectrumBe[i]>>distr[i];
    // cout<<i<<"\t"<<energyLoss[i] <<"\t"<<spectrum[i]<<endl;
  }
  simRead.close();

  // Cu data
  
  TGraph* ggCu = new TGraph(numberOfSimPoints, energyLoss, spectrumCu);

  ggCu->SetMarkerStyle(20);
  ggCu->SetMarkerSize(0.3);
  ggCu->SetMarkerColor(kGreen);
  ggCu->SetLineColor(kGreen);
  ggCu->SetLineStyle(1);
  ggCu->SetLineWidth(1);
  ggCu->SetFillColor(42);
  ggCu->SetTitle("Geant4 total");
  
  TGraphErrors* fnalCu = new TGraphErrors(7, pCuFNAL, sCuFNAL, xCuFNAL, yCuFNAL);
  fnalCu->SetMarkerStyle(20);
  fnalCu->SetMarkerColor(kGreen);
  fnalCu->SetLineColor(kGreen);
  // fnalCu->SetLineStyle(1);
  // fnalCu->SetLineWidth(0.001);
  fnalCu->SetFillColor(42);

  // U data
  
  TGraph* ggU = new TGraph(numberOfSimPoints, energyLoss, spectrumU);

  ggU->SetMarkerStyle(20);
  ggU->SetMarkerSize(0.3);
  ggU->SetMarkerColor(45);
  ggU->SetLineColor(45);
  ggU->SetLineStyle(1);
  ggU->SetLineWidth(1);
  ggU->SetFillColor(42);
  ggU->SetTitle("Geant4 total");
  
  TGraphErrors* fnalU = new TGraphErrors(7, pUFNAL, sUFNAL, xUFNAL, yUFNAL);
  fnalU->SetMarkerStyle(20);
  fnalU->SetMarkerColor(45);
  fnalU->SetLineColor(45);
  // fnalU->SetLineStyle(1);
  // fnalU->SetLineWidth(0.001);
  fnalU->SetFillColor(42);


  // Pb data
  
  TGraph* ggPb = new TGraph(numberOfSimPoints, energyLoss, spectrumPb);

  ggPb->SetMarkerStyle(20);
  ggPb->SetMarkerSize(0.3);
  ggPb->SetMarkerColor(41);
  ggPb->SetLineColor(41);
  ggPb->SetLineStyle(1);
  ggPb->SetLineWidth(1);
  ggPb->SetFillColor(42);
  ggPb->SetTitle("Geant4 total");
  
  TGraphErrors* fnalPb = new TGraphErrors(7, pPbFNAL, sPbFNAL, xPbFNAL, yPbFNAL);
  fnalPb->SetMarkerStyle(20);
  fnalPb->SetMarkerColor(41);
  fnalPb->SetLineColor(41);
  // fnalPb->SetLineStyle(1);
  // fnalPb->SetLineWidth(0.001);
  fnalPb->SetFillColor(42);

  // W data
  
  TGraph* ggW = new TGraph(numberOfSimPoints, energyLoss, spectrumW);

  ggW->SetMarkerStyle(20);
  ggW->SetMarkerSize(0.3);
  ggW->SetMarkerColor(6);
  ggW->SetLineColor(6);
  ggW->SetLineStyle(1);
  ggW->SetLineWidth(1);
  ggW->SetFillColor(42);
  ggW->SetTitle("Geant4 total");
  
  TGraphErrors* fnalW = new TGraphErrors(7, pWFNAL, sWFNAL, xWFNAL, yWFNAL);
  fnalW->SetMarkerStyle(20);
  fnalW->SetMarkerColor(6);
  fnalW->SetLineColor(6);
  // fnalW->SetLineStyle(1);
  // fnalW->SetLineWidth(0.001);
  fnalW->SetFillColor(42);

  // Fe data
  
  TGraph* ggFe = new TGraph(numberOfSimPoints, energyLoss, spectrumFe);

  ggFe->SetMarkerStyle(20);
  ggFe->SetMarkerSize(0.3);
  ggFe->SetMarkerColor(48);
  ggFe->SetLineColor(48);
  ggFe->SetLineStyle(1);
  ggFe->SetLineWidth(1);
  ggFe->SetFillColor(42);
  ggFe->SetTitle("Geant4 total");
  
  TGraphErrors* fnalFe = new TGraphErrors(7, pFeFNAL, sFeFNAL, xFeFNAL, yFeFNAL);
  fnalFe->SetMarkerStyle(20);
  fnalFe->SetMarkerColor(48);
  fnalFe->SetLineColor(48);
  // fnalFe->SetLineStyle(1);
  // fnalFe->SetLineWidth(0.001);
  fnalFe->SetFillColor(42);
  fnalFe->SetTitle("FNAL data");
  fnalFe->SetTitle("FNAL data");

  // C data
  
  TGraph* ggC = new TGraph(numberOfSimPoints, energyLoss, spectrumC);

  ggC->SetMarkerStyle(20);
  ggC->SetMarkerSize(0.3);
  ggC->SetMarkerColor(kRed);
  ggC->SetLineColor(kRed);
  ggC->SetLineStyle(1);
  ggC->SetLineWidth(1);
  ggC->SetFillColor(42);
  ggC->SetTitle("Geant4 total");
  
  TGraphErrors* fnalC = new TGraphErrors(7, pCFNAL, sCFNAL, xCFNAL, yCFNAL);
  fnalC->SetMarkerStyle(20);
  fnalC->SetMarkerColor(kRed);
  fnalC->SetLineColor(kRed);
  // fnalC->SetLineStyle(1);
  // fnalC->SetLineWidth(0.001);
  fnalC->SetFillColor(42);
  fnalC->SetTitle("FNAL data");

  // Be data

  TGraph* ggBe = new TGraph(numberOfSimPoints, energyLoss, spectrumBe);

  ggBe->SetMarkerStyle(20);
  ggBe->SetMarkerSize(0.3);
  ggBe->SetMarkerColor(kBlue);
  ggBe->SetLineColor(kBlue);
  ggBe->SetLineStyle(1);
  ggBe->SetLineWidth(1);
  ggBe->SetFillColor(42);
  ggBe->SetTitle("Geant4 total");
  
  TGraphErrors* fnalBe = new TGraphErrors(7, pBeFNAL, sBeFNAL, xBeFNAL, yBeFNAL);
  fnalBe->SetMarkerStyle(20);
  fnalBe->SetMarkerColor(kBlue);
  fnalBe->SetLineColor(kBlue);
  // fnalBe->SetLineStyle(1);
  // fnalBe->SetLineWidth(0.001);
  fnalBe->SetFillColor(42);
  fnalBe->SetTitle("FNAL data");

  // create a multigraph

  TMultiGraph* mg = 
  new TMultiGraph("mg", "n-Target total cross-sections");

  
  // add the graphs to the multigraph
  
  mg->Add(ggU); 
  mg->Add(fnalU); 
  
  mg->Add(ggPb); 
  mg->Add(fnalPb); 
  
  mg->Add(ggW); 
  mg->Add(fnalW); 
  
  mg->Add(ggCu); 
  mg->Add(fnalCu); 

  mg->Add(ggFe); 
  mg->Add(fnalFe); 

  mg->Add(ggC); 
  mg->Add(fnalC); 

  mg->Add(ggBe); 
  mg->Add(fnalBe); 

  TCanvas* myc = new TCanvas("myc", "n-Target total cross-sections");

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
  

 // draw the legend

  TLegend* legend = new TLegend(0.6,0.65,0.88,0.85);

  legend->SetTextFont(72);
  legend->SetTextSize(0.03);

  legend->AddEntry(fnalU,"U FNAL data","lpe1");
  legend->AddEntry(ggU,"G-G total U","l");

  legend->AddEntry(fnalPb,"Pb FNAL data","lpe1");
  legend->AddEntry(ggPb,"G-G total Pb","l");

  legend->AddEntry(fnalW,"W FNAL data","lpe1");
  legend->AddEntry(ggW,"G-G total W","l");

  legend->AddEntry(fnalCu,"Cu FNAL data","lpe1");
  legend->AddEntry(ggCu,"G-G total Cu","l");

  legend->AddEntry(fnalFe,"Fe FNAL data","lpe1");
  legend->AddEntry(ggFe,"G-G total Fe","l");

  legend->AddEntry(fnalC,"C FNAL data","lpe1");
  legend->AddEntry(ggC,"G-G total C","l");

  legend->AddEntry(fnalBe,"Be FNAL data","lpe1");
  legend->AddEntry(ggBe,"G-G total Be","l");
 
 

  legend->Draw();
  gr1->Draw("e1p");


}

