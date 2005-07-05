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
//    **********************************
//    *                                *
//    *      Tst51AnalysisManager.hh  *
//    *                                *
//    **********************************
// 

// the class Analysis creates and managed histograms and ntuples
///
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
 
#ifdef G4ANALYSIS_USE
#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
# include <AIDA/AIDA.h>

namespace AIDA 
{
  class ITree;
  class IHistogramFactory;
  class IAnalysisFactory;
  class IDataPoint;
  class ITupleFactory;
  class ITuple;
  class ITreeFactory;
};


class Tst51AnalysisManager { 

public:
  
  ~Tst51AnalysisManager();
  static Tst51AnalysisManager* getInstance();
  void book();
  void angularDistributionTransmittedGamma(G4double);
  void energyDistributionTransmittedGamma(G4double);
  void angularDistributionBackGamma(G4double);
  void energyDistributionBackGamma(G4double);
  void angularDistributionPostStep(G4double);
  void energyDistributionPostStep(G4double);
  void fillNtuple(G4double, G4double);
  void finish();

private:

  static Tst51AnalysisManager* instance;
  Tst51AnalysisManager();

  AIDA::IAnalysisFactory*  aFact; 
  AIDA::ITreeFactory*      treeFact;
  AIDA::ITree*             theTree;

  AIDA::IHistogramFactory*     histogramFactory;

  AIDA::IHistogram1D* h1;
  AIDA::IHistogram1D* h2;
  AIDA::IHistogram1D* h3;
  AIDA::IHistogram1D* h4;
  AIDA::IHistogram1D* angularPostStep;
  AIDA::IHistogram1D* energyPostStep;
  AIDA::ITupleFactory *tupFact;
  AIDA::ITuple *ntuple;

  //AIDA::ITupleFactory* tupFact;
  //  AIDA::ITuple *ntuple;
};
#endif
#endif




