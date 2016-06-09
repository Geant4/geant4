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
// Code developed by:
// S.Guatelli
//
//
// $Id: BrachyAnalysisManager.hh,v 1.3 2006/06/29 17:30:48 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
  void FillHistogramWithEnergy(G4double,G4double,G4double);
  void PrimaryParticleEnergySpectrum(G4double);
  void DoseDistribution(G4double,G4double);
  void finish();

private:
  static BrachyAnalysisManager* instance;

private:
  AIDA::IAnalysisFactory*  aFact;
  AIDA::ITree*             theTree;
  AIDA::IHistogramFactory *histFact;
  AIDA::ITreeFactory      *treeFact;
  AIDA::IHistogram2D *h1;
  AIDA::IHistogram1D *h2;
  AIDA::IHistogram1D *h3;
};
#endif
#endif



