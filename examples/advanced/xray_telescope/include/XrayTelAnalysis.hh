//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: XrayTelAnalysis.hh,v 1.6 2002-11-19 18:03:14 santin Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#ifdef G4ANALYSIS_USE_PLOTTER
  void plotAll();
#endif
  static XrayTelAnalysis* instance;

#ifdef G4ANALYSIS_USE
  AIDA::IAnalysisFactory  *analysisFactory;
  AIDA::ITree             *tree;
  AIDA::IHistogramFactory *histoFactory;
  AIDA::ITupleFactory     *tupleFactory;
#ifdef G4ANALYSIS_USE_PLOTTER
  AIDA::IPlotterFactory   *plotterFactory;
  AIDA::IPlotter          *plotter;
#endif
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

  //  G4std::ofstream asciiFile;

};

#endif 
