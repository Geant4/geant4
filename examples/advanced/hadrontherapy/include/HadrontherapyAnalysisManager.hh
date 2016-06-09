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
#ifdef G4ANALYSIS_USE
#ifndef HADRONTHERAPYANALYSISMANAGER_HH
#define HADRONTHERAPYANALYSISMANAGER_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
# include <AIDA/AIDA.h>

namespace AIDA{
  class ITree; 
  class IDataPoint;
  class IAnalysisFactory;
  class ITreeFactory;
};

class HadrontherapyAnalysisManager
{
private:
  HadrontherapyAnalysisManager();

public:
  ~HadrontherapyAnalysisManager();
  static HadrontherapyAnalysisManager* getInstance();
  
  void book();
  
  void Energy_Dep(G4double, G4double);
  void Energy_Event(G4int, G4double);
  
  void finish();

private:
  static HadrontherapyAnalysisManager* instance;

private:
  AIDA::IAnalysisFactory*  aFact;
  AIDA::ITree*             theTree;
  //  AIDA::ITreeFactory      *treeFact; 
  AIDA::IHistogramFactory *histFact;
  AIDA::ITupleFactory     *tupFact;
  AIDA::IHistogram1D *h1;
  AIDA::ITuple *ntuple;
};

#endif
#endif



