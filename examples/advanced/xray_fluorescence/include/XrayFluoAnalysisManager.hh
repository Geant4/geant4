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
// $Id: XrayFluoAnalysisManager.hh
// GEANT4 tag $Name: xray_fluo-V04-01-03 
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  30 Nov 2001 Guy Barrand : migrate to AIDA-2.2.
//  28 Nov 2001 Elena Guardincerri   Created
//  29 Nov 2002 mgration to AIDA 3.0 (Alfonso.mantero@ge.infn.it)
//
// -------------------------------------------------------------------

#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

//#ifndef XrayFluoAnalysisManager_h
//#define XrayFluoAnalysisManager_h 1

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"
#include "XrayFluoDataSet.hh"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "XrayFluoAnalysisMessenger.hh"

class G4Step;
//class XrayFluoEventAction;

namespace AIDA {
class IAnalysisFactory;
class IHistogramFactory;
class ITree;
class ITupleFactory;
class ITuple;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoAnalysisManager
{
public:
 
  virtual ~XrayFluoAnalysisManager();
  
  void book();
  
  void finish();
  
  //fill histograms with data from XrayFluoSteppingAction
  void analyseStepping(const G4Step* aStep);
  
 //fill histograms with data from XrayFluoEventAction
  void analyseEnergyDep(G4double eDep);
  
 //fill histograms with data from XrayFluoPrimarygeneratorAction
  void analysePrimaryGenerator(G4double energy);
  
  //method to call to create an instance of this class
  static XrayFluoAnalysisManager* getInstance();

  //method intended to chenge the name of the hbook output file
  void SetOutputFileName(G4String);

 
private:
  //private constructor in order to create a singleton
 

  XrayFluoAnalysisManager();
 
  G4String outputFileName;

  static XrayFluoAnalysisManager* instance;
  
  //pointer to the analysis messenger
  XrayFluoAnalysisMessenger* analisysMessenger;

  //XrayFluoEventAction* pEvent;

  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITree* tree;

  AIDA::IHistogramFactory *histogramFactory;
  //  ITupleFactory* tupleFactory;

  //  ITuple* tuple;

  AIDA::IHistogram1D*   histo_1;
  AIDA::IHistogram1D*   histo_2;
  AIDA::IHistogram1D*   histo_3;
  AIDA::IHistogram1D*   histo_4;
  AIDA::IHistogram1D*   histo_5;
  AIDA::IHistogram1D*   histo_6;
  AIDA::IHistogram1D*   histo_7;
  AIDA::IHistogram1D*   histo_8;
  AIDA::IHistogram1D*   histo_9;
  AIDA::IHistogram1D*   histo_10;
  AIDA::IHistogram1D*   histo_12;
  AIDA::IHistogram1D*   histo_11;
};

#endif



