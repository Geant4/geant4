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
#ifdef  G4ANALYSIS_USE
#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

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

};
#endif
#endif




