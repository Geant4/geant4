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
// Author: A. Pfeiffer (Andreas.Pfeiffer@cern.ch) 
//         (copy of his UserAnalyser class)
//
// History:
// -----------
//  7 Nov 2001   MGP  Implemented according to A. Pfeiffer's instructions
// 18 Nov 2001   G.Santin GammaRayTel analysis management modified
//               according to the new design
//
// -------------------------------------------------------------------
// Class description:
// Example of analysis in a simulation application (histograms, ntuples etc.)
// This class follows the singleton design pattern; 
// it is responsible for the analysis management and algorithms 
// Histograms are compliant with AIDA, except for the usage of IHistoManager, 
// which is Anaphe/Lizard-specific
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
#ifdef  G4ANALYSIS_USE

#ifndef GammaRayTelAnalysis_h
#define GammaRayTelAnalysis_h

#include "globals.hh"

// Histogramming from AIDA 
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IHistogram2D.h"

// Histogramming from Anaphe
#include "Interfaces/IHistoManager.h"

// Ntuples from Anaphe
#include "NtupleTag/LizardNTupleFactory.h"
#include "NtupleTag/LizardQuantity.h"

// Vectors from ?
#include "Interfaces/IVector.h"
#include "Interfaces/IVectorFactory.h"

// Plotting from Anaphe?
#include "Interfaces/IPlotter.h"
//#include "AIDA_Plotter/AIDAPlotter.h"

class G4Track;
class NTuple;
class GammaRayTelAnalysisMessenger;
class GammaRayTelDetectorConstruction;

class GammaRayTelAnalysis
{
public:

  ~GammaRayTelAnalysis();

  void BeginOfRun();
  void EndOfRun(G4int n);
  void EndOfEvent(G4int flag);

  void SetHisto1DDraw(G4String str) {histo1DDraw = str;};
  void SetHisto1DSave(G4String str) {histo1DSave = str;};
  void SetHisto2DDraw(G4String str) {histo2DDraw = str;};
  void SetHisto2DSave(G4String str) {histo2DSave = str;};
  void SetHisto2DMode(G4String str) {histo2DMode = str;};
  G4String GetHisto2DMode() {return histo2DMode;};

  void InsertPositionXZ(double x, double z);
  void InsertPositionYZ(double y, double z);
  void InsertEnergy(double en);
  void InsertHits(int nplane);

  static GammaRayTelAnalysis* getInstance();

private:

  GammaRayTelAnalysis();

  void plot1D(IHistogram1D* histo);
  void plot2D(IHistogram2D* histo);
  void Plot();

  static GammaRayTelAnalysis* instance;

  IHistoManager* histoManager;
  IVectorFactory* vectorFactory;
  IPlotter* plotter;

  // ---- NOTE ----
  // Histograms are compliant to AIDA interfaces, ntuples are Lizard specific
 
  Lizard::NTuple* ntuple;
  Lizard::NTupleFactory* ntFactory;

  // Quantities for the ntuple
  Lizard::Quantity<float> ntEnergy;
  Lizard::Quantity<int>   ntPlane;
  Lizard::Quantity<float> ntX;
  Lizard::Quantity<float> ntY;
  Lizard::Quantity<float> ntZ;

  GammaRayTelDetectorConstruction*    GammaRayTelDetector;

  G4String histo1DDraw;
  G4String histo1DSave;
  G4String histo2DDraw;
  G4String histo2DSave;
  G4String histo2DMode;

  GammaRayTelAnalysisMessenger* analysisMessenger;

};

#endif 
#endif 
