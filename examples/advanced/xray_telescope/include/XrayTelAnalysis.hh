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
//
// $Id$
//
// Author: A. Pfeiffer (Andreas.Pfeiffer@cern.ch) 
//         (copy of his UserAnalyser class)
//
// History:
// -----------
//  7 Nov 2001   MGP  Implemented according to A. Pfeiffer's instructions
//
// -------------------------------------------------------------------
// Class description:
// Example of analysis in a simulation application (histograms, ntuples etc.)
// This class follows the singleton design pattern; 
// it is responsible for the analysis management and algorithms 
// Histograms are compliant with AIDA 1.0 (the header files are
// copies from the Anaphe Interfaces directory), except for the usage of
// IHistoManager which is Anaphe/Lizard-specific
// For ntuples, for which the AIDA interfaces and compliant implementations
// are still in progress at the time of the current Geant4 release, 
// an implementation with Anaphe/Lizard is shown
// Other implementations specific to an analysis system are possible too
// (see, for instance, JAS and OpenScientist documentation from the links
// in http://aida.freehep.org/)
// The implementation of the usage of ntuples shown in this example
// is expected to change in future Geant4 releases, when AIDA interfaces 
// and related implementations would be available
// Further documentation is available from: http://www.ge.infn.it/geant4/lowE/
//                                          http://aida.freehep.org/
//                                          http://cern.ch/anaphe/

// -------------------------------------------------------------------

#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include "G4ios.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/IAnalysisFactory.h"

#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"

#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IHistogram3D.h"

#include "AIDA/IPlotterFactory.h"
#include "AIDA/IPlotterRegion.h"
#include "AIDA/IPlotter.h"

#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

#include "AIDA/IManagedObject.h"

// Histogramming from AIDA 
class IAnalysisFactory;
class ITree;
class IHistogramFactory;
class ITupleFactory;
class ITuple;
#endif

class G4Track;

class XrayTelAnalysis
{
public:

  ~XrayTelAnalysis();

  void book();
  
  void finish();
  
  void analyseStepping(const G4Track& track, G4bool entering);

  static XrayTelAnalysis* getInstance();

private:

  XrayTelAnalysis();

// #ifdef G4ANALYSIS_USE_PLOTTER
//   void plotAll();
// #endif
  static XrayTelAnalysis* instance;

#ifdef G4ANALYSIS_USE
  AIDA::IAnalysisFactory  *analysisFactory;
  AIDA::ITree             *tree;
  AIDA::IHistogramFactory *histoFactory;
  AIDA::ITupleFactory     *tupleFactory;
// #ifdef G4ANALYSIS_USE_PLOTTER
//   AIDA::IPlotterFactory   *plotterFactory;
//   AIDA::IPlotter          *plotter;
// #endif
#endif
  // Quantities for the ntuple
  G4double eKin;
  G4double x;
  G4double y;
  G4double z;
  G4double dirX;
  G4double dirY;
  G4double dirZ;

  G4String asciiFileName;
  G4String histFileName;
  G4String histFileType;

#ifdef G4ANALYSIS_USE
  AIDA::IHistogram1D *h1;
  AIDA::IHistogram2D *h2;
  AIDA::IHistogram1D *h3;
  AIDA::IHistogram2D *h4;
  AIDA::ITuple * ntuple;
#endif
  //  std::ofstream asciiFile;

};

#endif 
