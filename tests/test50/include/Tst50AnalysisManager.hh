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
// $Id: Tst50AnalysisManager.hh,v 1.6 2003-01-28 08:57:49 guatelli Exp $
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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class Tst50PrimaryGeneratorAction;
class Tst50AnalysisManager
{
public:

~Tst50AnalysisManager();

void book();

void finish();

static Tst50AnalysisManager* getInstance();

void energy_deposit(G4double);
  void energy_transmitted(G4double);
  void fill_data(G4double,G4double, G4double);
  void energy_backscatter(G4double);
  void angleB(G4double );
 void angleT(G4double );

private:
static Tst50AnalysisManager* instance;

private:
Tst50AnalysisManager();

private:

AIDA::IAnalysisFactory*  aFact;
AIDA::ITree*             theTree;
AIDA::IHistogramFactory *histFact;
AIDA::ITupleFactory     *tupFact;
AIDA::ITreeFactory      *treeFact;
AIDA::ITuple *ntuple;
AIDA::IHistogram1D *h1;
AIDA::IHistogram1D *h2;
AIDA::IHistogram1D *h3;
AIDA::IHistogram1D *h4;
AIDA::IHistogram1D *h5;
private:
  G4double initial_energy;
  Tst50PrimaryGeneratorAction* p_Primary;
  G4double en;
  G4double range;
};
#endif
#endif




