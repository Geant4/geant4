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
// $Id: Tst50AnalysisManager.hh,v 1.14 2003-05-17 11:59:35 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    **********************************
//    *                                *
//    *      BrachyAnalysisManager.hh  *
//    *                                *
//    **********************************
// 

// the class Analysis creates and managed histograms and ntuples
//
#ifdef G4ANALYSIS_USE
#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include "g4std/vector"
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


class Tst50AnalysisManager { 

public:
  
  ~Tst50AnalysisManager();
  static Tst50AnalysisManager* getInstance();
  void book();
  void attenuation_coeffiecient(G4int,G4double,G4double,G4double);
  void StoppingPower(G4int,G4double,G4double);
  void CSDARange(G4int,G4double,G4double);
  void trasmission(G4int,G4double,G4double,G4double,G4double,G4double); 
  void finish();

private:

  static Tst50AnalysisManager* instance;
  Tst50AnalysisManager();

  AIDA::IAnalysisFactory*  aFact; 
  AIDA::ITreeFactory*      treeFact;
  AIDA::ITree*             theTree;
  AIDA::IDataPointSetFactory * dpsf;
  AIDA::IDataPointSet * dpsa;  
  AIDA::IDataPointSet * dpsa1;
  AIDA::IDataPointSet * dpsa2;
  AIDA::IDataPointSet * dpsa3;
};
#endif
#endif




