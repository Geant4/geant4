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
};
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

public:

  AIDA::IAnalysisFactory*  aFact;
  AIDA::ITree*             theTree;
  AIDA::IHistogramFactory *histFact;
  AIDA::ITreeFactory      *treeFact;
  AIDA::IHistogram1D *h1;
  AIDA::IHistogram1D *h2;
};

#endif

#endif
