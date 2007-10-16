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
// $Id: TestDens.cc,v 1.3 2007-10-16 20:37:44 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestDens.hh"

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

ClassImp(TestDens)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test dens</h1>

     <p>
     This class tests method G4Incl::dens
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestDens.png"/>
     </p>
     END_HTML */
  
TestDens::TestDens()
{
  setUnitName("G4Incl::Dens");
  setOriginalUnitName("dens");
  setPlotFileName("htmldoc/plots/TestDens.png");
  setLogFileName("htmldoc/logs/TestDens.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestDens::~TestDens()
{

}

void TestDens::runMe()
{
  // This handles the actual testing.

  G4Incl *incl = new G4Incl();
  
  const Float_t errorMarginal = 1e-6;

  Float_t r0 = 0.1;
  Float_t adif = 1.0;
  Float_t integral = 0.0;
  
  // Set up arrays for data collection:
  const Int_t points = 100;
  const Float_t r_step = 0.05;
  Float_t r[points];
  r[0] = 0.0;
  Float_t cpp_f_r[points];
  Float_t fort_f_r[points];
  Float_t relativeDifference[points];

  // Data for struct dton (common gDton):
  // Values are from code module init_incl42.cc
  const int dtonsize = 13;
  float cData[dtonsize] = {0.88688076e+00,-0.34717093e+00,-.30502380e+01,
			   .56207766e+02,-.74957334e+03,.53365279e+04,-.22706863e+05,
			   .60434469e+05,-.10292058e+06,.11223357e+06,-.75925226e+05,
			   .29059715e+05,-.48157368e+04};

  float dData[dtonsize] = {.23135193e-01,-.85604572e+00,.56068193e+01,
			   -.69462922e+02,.41631118e+03,-.12546621e+04,.12387830e+04,
			   .33739172e+04,-.13041151e+05,.19512524e+05,-.15634324e+05,
			   .66231089e+04,-.11698185e+04};

  // C++:
  G4Dton *dton = (G4Dton*) malloc(sizeof(G4Dton));
  G4Ws *ws = (G4Ws*) malloc(sizeof(G4Ws));
  ws->r0 = r0;
  ws->adif = adif;

  // FORTRAN:
  gWs->r0 = r0;
  gWs->adif = adif;

  //Set initial values to common gDton and struct dton:
  for(Int_t i = 0; i < dtonsize; i++) {
    dton->c[i] = cData[i];
    gDton->c[i] = cData[i];
		
    dton->d[i] = dData[i];
    gDton->d[i] = dData[i];

    dton->fn = 1.0;
    gDton->fn = 1.0;
  }

  incl->setDtonData(dton);
  incl->setWsData(ws);
	
  for(Int_t i = 1; i < points; i++) {
    r[i] = r[i-1] + r_step;	
    cpp_f_r[i] = incl->dens(r[i]);
    fort_f_r[i] = dens_(&r[i]);
    relativeDifference[i] = 100.0*(cpp_f_r[i] - fort_f_r[i])/fort_f_r[i];
    integral = integral + abs(cpp_f_r[i] - fort_f_r[i]);
  }

  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(points, r, cpp_f_r);
  TGraph *fort_graph = new TGraph(points, r, fort_f_r);

  fort_graph->GetXaxis()->SetTitle("q");
  fort_graph->GetYaxis()->SetTitle("dens(q)");
  fort_graph->SetTitle("Function dens");

  fort_graph->Draw("al");
  cpp_graph->Draw("p, same");	

  c1->cd(4);
  TGraph *diff_graph = new TGraph(points, r, relativeDifference);
  diff_graph->SetTitle("Relative difference C++/FORTRAN");
  diff_graph->GetXaxis()->SetTitle("q");
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
    setTestStatus(true);
  }
  else {
    setTestStatus(false);
  }
}
