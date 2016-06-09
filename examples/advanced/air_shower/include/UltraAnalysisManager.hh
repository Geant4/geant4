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
//
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraAnalysisManager.hh
//    ****************************************************
//
//    Class used for analysis procedures if the environment variable 
//    G4ANALYSIS_USE is set. AIDA is supported.
//
#ifdef G4ANALYSIS_USE 
#ifndef ULTRAANALYSISMANAGER_HH
#define ULTRAANALYSISMANAGER_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

#ifdef G4ANALYSIS_USE
namespace AIDA{
  class ITree;
  class IHistogramFactory;
  class IAnalysisFactory;
  class ITreeFactory;
}
#endif

class UltraAnalysisManager
{
private:
  UltraAnalysisManager();

public:

  ~UltraAnalysisManager();
  static UltraAnalysisManager* getInstance();
  void book();
  void FillHistogram(G4int i, G4double f);
  void finish();

  
private:
  static UltraAnalysisManager* instance;

  AIDA::IAnalysisFactory*  aFact;
  AIDA::ITree*             theTree;
  AIDA::IHistogramFactory *histFact;
  AIDA::ITreeFactory      *treeFact;
  AIDA::IHistogram1D *h1;
  AIDA::IHistogram1D *h2;
};

#endif

#endif
