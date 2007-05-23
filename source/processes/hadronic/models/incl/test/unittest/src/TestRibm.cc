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
// $Id: TestRibm.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestRibm.hh"

#include "InclAblaTestDifferencePlotter.hh"

#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"

//#include "InclAblaFunctions.hh"
#include "functionwrapper.hh"
#include "commonwrapper.hh"

#ifndef __CINT__
#include "G4Incl.hh"
#include "G4InclDataDefs.hh"
#endif

#include <iostream>

ClassImp(TestRibm)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Random number generator test</h1>

     <p>
     This is a test for the basic random number generator included in INCL4.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestRibm.png"/>
     </p>
     END_HTML */
  
TestRibm::TestRibm()
{
  setUnitName("G4Incl::standardRandom");
  setOriginalUnitName("ribm");
  setPlotFileName("htmldoc/plots/TestRibm.png");
  setLogFileName("htmldoc/logs/TestRibm.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestRibm::~TestRibm()
{

}

void TestRibm::runMe()
{
  // This handles the actual testing.
	
  std::cout <<"Testing Ribm...." << std::endl;

  G4Incl *incl = new G4Incl();
  
  const Int_t collectNumbers = 100000;

  // Set the error marginal:
  const Float_t errorMarginal = 1e-9;
	  
  // Histogramming settings:
  const Int_t histogramBins = 100;
		          
  // Set same random seeds for both C++ and FORTRAN random number generators.
  int seed = 38457;
  double rndm; //Random numbers

  // For FORTRAN code
  int fseed = seed;
  float rndmf;

  // Histograms for random numbers:
  // C++ routine:
  TH1D *cpprndm = new TH1D("cpprndm", "Subroutine ribm", histogramBins, 0.0, 1.0);
  // FORTRAN routine
  TH1D *frndm = new TH1D("frndm", "Subroutine ribm", histogramBins, 0.0, 1.0);  
  // Relative difference:
  TH1D *relativeDiff = new TH1D("relativeDiff", "Relative differences between C++ and FORTRAN code", histogramBins, 0.0, 1.0);

  // Collect random numbers generated both codes and store them into histograms:
  for(Int_t i = 0; i < collectNumbers; i++) {
    incl->standardRandom(&rndm, &seed);
    ribm_(&rndmf, &fseed);
    cpprndm->Fill(rndm);
    frndm->Fill(rndmf);
  }

  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);

  cpprndm->SetMinimum(0);
  cpprndm->SetMaximum(1500);
  cpprndm->GetXaxis()->SetTitle("Random number");

  // Draw plot showing both histograms:
  cpprndm->Draw();
  frndm->Draw("P,same");

  c1->cd(4);
  //  TGraph *diff = InclAblaTestDifferencePlotter::RelativeDiff(cpprndm, frndm, "Random number");
  //  diff->Draw("ap");

    // Calculate the difference of the histograms.
  Int_t integral = 0;
  TGraph *relDiff = new TGraph(cpprndm->GetNbinsX());
  for(Int_t i = 0; i < cpprndm->GetNbinsX(); i++) {
    if(frndm->GetBinContent(i) != 0) {
      relativeDiff->SetBinContent(i, (100.0*(cpprndm->GetBinContent(i) - frndm->GetBinContent(i))/frndm->GetBinContent(i)));
    }
    else {
      relativeDiff->SetBinContent(i, 0);
    }
    integral = integral + abs(cpprndm->GetBinContent(i) - frndm->GetBinContent(i));
    if(frndm->GetBinContent(i) != 0) {
      relDiff->SetPoint(i, cpprndm->GetBinLowEdge(i), 100.0*(cpprndm->GetBinContent(i) - frndm->GetBinContent(i))/frndm->GetBinContent(i));
    }
    else {
      relDiff->SetPoint(i, cpprndm->GetBinLowEdge(i), 0.0);
    }
  }
  
  relDiff->SetTitle("Relative difference C++/FORTRAN");
  relDiff->GetXaxis()->SetTitle("Random number");
  relDiff->GetYaxis()->SetTitle("Relative difference (%)");
  relDiff->Draw("ap");

  c1->SaveAs(getPlotFileName());

  std::cout <<"Integral: " << integral << std::endl;

  delete c1;
  delete cpprndm;
  delete frndm;
  delete relativeDiff;

  // The integral over the difference of the bins of the histograms should be 0
  // if both C++ and FORTRAN results are the same.
  if(integral < errorMarginal) {
    // Passed the test
    setTestStatus(true);
  }
  else {
    // Failed the test
    setTestStatus(false);
  }
}
