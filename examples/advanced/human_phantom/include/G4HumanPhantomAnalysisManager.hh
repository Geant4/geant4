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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//

#ifdef  G4ANALYSIS_USE
#ifndef G4HUMANPHANTOMANALYSISMANAGER_HH
#define G4HUMANPHANTOMANALYSISMANAGER_HH 1

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include <AIDA/AIDA.h>

namespace AIDA 
{
  class ITree;
  class IHistogramFactory;
  class IAnalysisFactory;
  class ITreeFactory;
};

class G4HumanPhantomAnalysisManager { 

public:
  
  ~G4HumanPhantomAnalysisManager();
  static G4HumanPhantomAnalysisManager* getInstance();

  void book();
  void particlePath(G4double);
  void particleProjectionXY(G4double, G4double);
  void particleProjectionYZ(G4double, G4double);
  void particleProjectionZX(G4double, G4double);
  void bodypartEnergyDep(G4double,G4double);
  void innerBreastEnergyDep(G4int, G4int, G4double);
  void finish();

private:
  static G4HumanPhantomAnalysisManager* instance;
  G4HumanPhantomAnalysisManager();

  AIDA::IAnalysisFactory*  aFact; 
  AIDA::ITreeFactory*      treeFact;
  AIDA::ITree*             theTree;
  AIDA::IHistogramFactory* histogramFactory;
  AIDA::IHistogram1D*      histogramParticlePath;
  AIDA::IHistogram2D*      projectionXY;
  AIDA::IHistogram2D*      projectionYZ;
  AIDA::IHistogram2D*      projectionZX;
  AIDA::IHistogram2D*      energy;
  AIDA::IHistogram2D*      innerBreast;
};
#endif
#endif




