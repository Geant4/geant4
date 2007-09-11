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
// $Id: TestSigReac.cc,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestSigReac.hh"

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

#include "Riostream.h"

ClassImp(TestSigReac)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test sig_reac</h1>

     <p>
     This class tests method G4Incl::crossSection.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestSigReac.png"/>
     </p>
     END_HTML */
  
TestSigReac::TestSigReac()
{
  setUnitName("G4Incl::crossSection");
  setOriginalUnitName("sig_reac");
  setPlotFileName("htmldoc/plots/TestSigReac.png");
  setLogFileName("htmldoc/logs/TestSigReac.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestSigReac::~TestSigReac()
{

}

void TestSigReac::runMe()
{
  // This handles the actual testing.

  const double errorMarginal = 1e-3;
  
  G4Incl *incl = new G4Incl();

  double E = 50.0;
  double minA = 28.0;
  const int nuclei = 50;

  double A = minA;

  double Ai[nuclei];
  double cppCrossSection[nuclei];
  double fortCrossSection[nuclei];
  float fortCX;
  int bullet = 1; // Proton

  double relativeDifference[nuclei];
  double integral = 0.0;
  
  for(Int_t i = 0; i < nuclei; i++) {
    Ai[i] = A;
    cppCrossSection[i] = incl->crossSection(bullet, E, Ai[i]);
    sig_reac__(&bullet, &E, &Ai[i], &fortCX);
    fortCrossSection[i] = fortCX;
    //    cout << cppCrossSection[i] << " \t " << fortCX << endl;
    relativeDifference[i] = 100.0*(cppCrossSection[i] - fortCrossSection[i])/fortCrossSection[i];
    integral = integral + abs(cppCrossSection[i] - fortCrossSection[i]);

    A = A + 1.0;
  }
  
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(nuclei, Ai, cppCrossSection);
  TGraph *fort_graph = new TGraph(nuclei, Ai, fortCrossSection);

  fort_graph->GetXaxis()->SetTitle("A");
  fort_graph->GetYaxis()->SetTitle("sig_reac(bullet=proton, E=50MeV, A)");
  fort_graph->SetTitle("Function sig_reac");

  fort_graph->Draw("al");
  cpp_graph->Draw("p, same");	

  c1->cd(4);
  TGraph *diff_graph = new TGraph(nuclei, Ai, relativeDifference);
  diff_graph->SetTitle("Relative difference C++/FORTRAN");
  diff_graph->GetXaxis()->SetTitle("r");
  diff_graph->GetYaxis()->SetTitle("Relative difference (%)");

  diff_graph->Draw("ap");

  c1->SaveAs(getPlotFileName());

  // Clean up...
  delete c1;
  delete cpp_graph;
  delete fort_graph;
  delete diff_graph;

  // The integral over the difference of the bins of the histograms should be
  // smaller than errorMarginal if both C++ and FORTRAN results are the same.
  if(integral < errorMarginal) {
    // Passed the test
    setTestStatus(true);
  }
  else {
    // Failed the test
    setTestStatus(false);
  }
}
