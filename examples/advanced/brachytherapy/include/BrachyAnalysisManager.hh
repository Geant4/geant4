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
// $Id: BrachyAnalysisManager.hh,v 1.5 2002-11-18 15:18:35 guatelli Exp $
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

class BrachyAnalysisManager
{
public:

~BrachyAnalysisManager();

void book();

void finish();

static BrachyAnalysisManager* getInstance();

void fill_Tuple(G4double,G4double,G4double,G4float);
void hist(G4double,G4double,G4float);
void Spectrum(G4double);

private:

G4double xx,zz,yy;
G4float  en; 
G4double  x,y,z;
static BrachyAnalysisManager* instance;

private:
BrachyAnalysisManager();

private:

AIDA::IAnalysisFactory*  aFact;
AIDA::ITree*             theTree;
AIDA::IHistogramFactory *histFact;
AIDA::ITupleFactory     *tupFact;
AIDA::ITreeFactory      *treeFact;
AIDA::IHistogram2D *h1;
AIDA::IHistogram1D *h2;
AIDA::ITuple *ntuple;
};

#endif




