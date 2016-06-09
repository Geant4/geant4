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
// $Id: MedLinacAnalysisManager.hh,v 1.3 2004/05/14 18:25:39 mpiergen Exp $
//
//
// Code developed by: M. Piergentili
//
// the class Analysis creates and managed histograms
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
  class ITreeFactory;
};

class MedLinacAnalysisManager
{
private:
  MedLinacAnalysisManager();

public:

  ~MedLinacAnalysisManager();
  static MedLinacAnalysisManager* getInstance();
  void book();
  void FillHistogram1WithEnergy(G4double,G4double,G4float);
  void FillHistogram3WithEnergy(G4double,G4double,G4float);
  void FillHistogram4WithEnergy(G4double,G4float);
  void FillHistogram5WithEnergy(G4double,G4float);
  void FillHistogram6WithEnergy(G4double,G4float);
  void FillHistogram7WithEnergy(G4double,G4float);
  void FillHistogram8WithEnergy(G4double,G4float);
  void PrimaryParticleEnergySpectrum(G4double);
  void finish();

private:

  G4double xx,zz,yy;
  G4float  en; 
  G4double  x,y,z;
  static MedLinacAnalysisManager* instance;

private:

  AIDA::IAnalysisFactory*  aFact;
  AIDA::ITree*             theTree;
  AIDA::IHistogramFactory *histFact;
  AIDA::ITreeFactory      *treeFact;
  AIDA::IHistogram2D *h1;
  AIDA::IHistogram1D *h2;
  AIDA::IHistogram2D *h3;
  AIDA::IHistogram1D *h4;
  AIDA::IHistogram1D *h5;
  AIDA::IHistogram1D *h6;
  AIDA::IHistogram1D *h7;
  AIDA::IHistogram1D *h8;
};

#endif
#endif



