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
// Code developed by:
// S.Guatelli
//
//
// $Id: BrachyAnalysisManager.hh,v 1.9 2004/03/11 15:38:42 guatelli Exp $
// GEANT4 tag $Name: geant4-06-01 $
//
//    **********************************
//    *                                *
//    *      BrachyAnalysisManager.hh  *
//    *                                *
//    **********************************
// 
//
// the class Analysis creates and managed histograms and ntuples
//
#ifdef G4ANALYSIS_USE
#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IAnalysisFactory.h"

namespace AIDA{
  class ITree;
  class IHistogramFactory;
  class IAnalysisFactory;
  class ITupleFactory;
  class ITuple;
  class ITreeFactory;
};

class BrachyAnalysisManager
{
private:
  BrachyAnalysisManager();

public:

  ~BrachyAnalysisManager();
  static BrachyAnalysisManager* getInstance();
  void book();
  void FillNtupleWithEnergy(G4double,G4double,G4double,G4float);
  void FillHistogramWithEnergy(G4double,G4double,G4double);
  void PrimaryParticleEnergySpectrum(G4double);
  void DoseDistribution(G4double,G4double);
  void finish();

  
private:

  //  G4double xx,zz,yy;
  //G4float  en; 
  //G4double  x,y,z;
  static BrachyAnalysisManager* instance;

private:

  AIDA::IAnalysisFactory*  aFact;
  AIDA::ITree*             theTree;
  AIDA::IHistogramFactory *histFact;
  AIDA::ITupleFactory     *tupFact;
  AIDA::ITreeFactory      *treeFact;
  AIDA::IHistogram2D *h1;
  AIDA::IHistogram1D *h2;
  AIDA::IHistogram1D *h3;
  AIDA::ITuple *ntuple;
};

#endif
#endif



