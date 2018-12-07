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
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//     
//
//      ------------ GammaRayTelAnalysis  ------
//           by R.Giannitrapani, F. Longo & G.Santin (30 nov 2000)
//
// 03.04.2013 F.Longo and L.Pandola
// - migrated to G4AnalysisManager
//
// 07.12.2001 A.Pfeiffer
// - integrated Guy's addition of the ntuple
//
// 06.12.2001 A.Pfeiffer
// - updating to new design (singleton)
//
// 22.11.2001 G.Barrand
// - Adaptation to AIDA
// -------------------------------------------------------------------
// Class description:
// Example of analysis in a simulation application (histograms, ntuples etc.)
// This class follows the singleton design pattern; 
// it is responsible for the analysis management and algorithms 
//
// -------------------------------------------------------------------//

#ifndef GammaRayTelAnalysis_h 
#define GammaRayTelAnalysis_h 1

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

// uncomment g4xml.hh and comment g4root.hh for a XML-based output file

#include "g4root.hh"
//#include "g4xml.hh"


class GammaRayTelAnalysisMessenger;
class GammaRayTelDetectorConstruction;

class GammaRayTelAnalysis {
public:
  virtual ~GammaRayTelAnalysis();
  
public:


  //  void BeginOfRun(G4int n);

  void BeginOfRun();
  void EndOfRun();
  void EndOfEvent(G4int flag);

  void Init();
  void Finish();

  void SetHisto2DMode(G4String str) {histo2DMode = str;};
  G4String GetHisto2DMode() {return histo2DMode;};
  
  void InsertPositionXZ(double x, double z);
  void InsertPositionYZ(double y, double z);
  void InsertEnergy(double en);
  void InsertHits(int nplane);

  void setNtuple(float E, float p, float x, float y, float z);

  static GammaRayTelAnalysis* getInstance();

private:

  GammaRayTelAnalysis();

  //void plot1D(IHistogram1D* histo);
  //void plot2D(IHistogram2D* histo);
  void Plot();

private:
  static GammaRayTelAnalysis* instance;

  const GammaRayTelDetectorConstruction*    GammaRayTelDetector;

  G4String histo2DMode;
  G4String histoFileName;

  GammaRayTelAnalysisMessenger* analysisMessenger;
};


#endif
