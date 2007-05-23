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
// $Id: InclAblaTestDifferencePlotter.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "InclAblaTestDifferencePlotter.hh"

#include <iostream>

ClassImp(InclAblaTestDifferencePlotter)

//////////////////////////////////////////////////////////////////////////
/* BEGIN_HTML
<h1>InclAblaTestDifferencePlotter</h1>

<h2>Overview</h2>


<h2>Examples of usage</h2>

END_HTML */

InclAblaTestDifferencePlotter::InclAblaTestDifferencePlotter()
{

}

InclAblaTestDifferencePlotter::~InclAblaTestDifferencePlotter()
{

}

TGraph* InclAblaTestDifferencePlotter::RelativeDiff(TH1D *hist1, TH1D *hist2, TString XAxisLabel)
{

  return 0;
}

TGraph* InclAblaTestDifferencePlotter::RelativeDiff(TH1F *hist1, TH1F *hist2, TString XAxisLabel)
{
  // Compares two histograms, calculates their relative differences and returns the result as a graph.
  
  Int_t hist1Bins = hist1->GetNbinsX();
  Int_t hist2Bins = hist2->GetNbinsX();

  if(hist1Bins != hist2Bins) {
    std::cout <<"InclAblaTestDifferencePlotter error: Histograms must have equal number of bins! " << std::endl;
    return 0;
  }

  TGraph *differenceGraph = new TGraph(hist1Bins);

  for(Int_t i = 0; i < hist1Bins; i++) {
    differenceGraph->SetPoint(i, i, 100*((hist1->GetBinContent(i) - hist2->GetBinContent(i))/hist2->GetBinContent(i)));
  }

  differenceGraph->GetXaxis()->SetTitle(XAxisLabel);
  
  return (TGraph*) differenceGraph->Clone("differences");
}

TGraph* InclAblaTestDifferencePlotter::RelativeDiff(TH1D *hist1, TH1F *hist2, TString XAxisLabel)
{

  return 0;
}

TGraph* InclAblaTestDifferencePlotter::RelativeDiff(TGraph *graph1, TGraph *graph2, Int_t points, TString XAxisLabel)
{
  // Compares two graphs and calculates their relative difference.
  
  TGraph *differenceGraph = new TGraph(points);
  double x1, y1;
  double x2, y2;
  
  for(Int_t i = 0; i < points; i++) {
    graph1->GetPoint(i, x1, y1);
    graph2->GetPoint(i, x2, y2);
    differenceGraph->SetPoint(i, i, 100*((y1 - y2)/y2));
  }

  differenceGraph->SetTitle("Relative differences C++/FORTRAN");
  differenceGraph->GetXaxis()->SetTitle(XAxisLabel);
  differenceGraph->GetYaxis()->SetTitle("Relative difference (%)");

  return ((TGraph*) differenceGraph->Clone("diffGraph"));
}
