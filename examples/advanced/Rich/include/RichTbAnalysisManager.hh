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
// Rich advanced example for Geant4
// RichTbAnalysisManager.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbAnalysisManager_h
#define RichTbAnalysisManager_h 1

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"



#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif
#include "RichTbAnalysisMessenger.hh"

class G4Step;

#ifdef G4ANALYSIS_USE
namespace AIDA{
  class IAnalysisFactory;
  class IHistogramFactory;
  class ITree;
  class ITupleFactory;
  class ITuple;
}
#endif

class G4Run;
class G4Event;
class G4Step;
class G4Timer;

class G4SteppingManager;
class RichTbRunConfig;

class RichTbAnalysisManager{
public: 

  RichTbAnalysisManager();
  virtual ~RichTbAnalysisManager();

#ifdef G4ANALYSIS_USE

  void book();
  void finish();

  static RichTbAnalysisManager* getInstance();
  void SetOutputFileName(G4String);





#endif



  void BeginOfRunAnalysis(); 
  void EndOfRunAnalysis(); 
  void BeginOfEventAnalysis(const G4Event*); 
  void EndOfEventAnalysis(const G4Event*); 
  void StepAnalysis(const G4SteppingManager* );
    
  void bumpNumPhotonsBeforeAerogel() { NumPhotBeforeAerogel++ ; }
  G4int getNumPhotBeforeAerogel() {return  NumPhotBeforeAerogel; }

  void bumpNumPhotonsAfterFilter(){  NumPhotAfterFilter++ ; }
  G4int getNumPhotAfterFilter() {return  NumPhotAfterFilter ; }
  void bumpNumPhotonsBeforeFilter(){  NumPhotBeforeFilter++ ; }
  G4int getNumPhotBeforeFilter() {return  NumPhotBeforeFilter ; }
  void bumpNumPhotonsBeforeMirror(){  NumPhotBeforeMirror++ ; }
  G4int getNumPhotBeforeMirror() {return  NumPhotBeforeMirror ; }

  void bumpNumPhotonsAfterMirror(){  NumPhotAfterMirror++ ; }
  G4int getNumPhotAfterMirror() {return  NumPhotAfterMirror ; }

  void bumpNumPhotonsAtHpd1Input(){  NumPhotAtHpd1Input++ ; }
  G4int getNumPhotAtHpd1Input() {return  NumPhotAtHpd1Input ; }
  void bumpNumPhotonsAtHpd2Input(){  NumPhotAtHpd2Input++ ; }
  G4int getNumPhotAtHpd2Input() {return  NumPhotAtHpd2Input ; }
  void bumpNumPhotonsAtHpd3Input(){  NumPhotAtHpd3Input++ ; }
  G4int getNumPhotAtHpd3Input() {return  NumPhotAtHpd3Input ; }
  void bumpNumPhotonsAtHpd4Input(){  NumPhotAtHpd4Input++ ; }
  G4int getNumPhotAtHpd4Input() {return  NumPhotAtHpd4Input ; }
  void bumpNumHitInSi(){  NumHitInSi++ ; }
  G4int getNumHitInSi() {return  NumHitInSi ; }
  G4double getMeanHXCoord() {return  MeanHXCoord ; }
  G4double getMeanHYCoord() {return  MeanHYCoord ; }
  G4double getNumHXCoord() {return  NumHXCoord ; }
  G4double getNumHYCoord() {return  NumHYCoord ; }
  void setMeanHXCoord(G4double cx );
  void setMeanHYCoord(G4double cy );

public:

  G4String outputFileName;

  static RichTbAnalysisManager* instance;

#ifdef G4ANALYSIS_USE

  RichTbAnalysisMessenger* analisysMessenger;

  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITree* tree;
  AIDA::IHistogramFactory* histogramFactory;  

  AIDA::IHistogram1D* getfhistoNPhotG() {return fhistoNrPhotG;}  
  AIDA::IHistogram1D* getfhistoNBeforeMirror() {return fhistoNBeforeMirror;}  
  AIDA::IHistogram1D* getfhistoWBeforeMirror() {return fhistoWBeforeMirror;}  
  AIDA::IHistogram1D* getfhistoWAfterMirror() {return fhistoWAfterMirror;}

  AIDA::IHistogram1D* getfhistoCkvProdSmall() {return fhistoCkvProdSmall;}  
  AIDA::IHistogram1D* getfhistoEmisZ() {return fhistoEmisZ;}  
 
  AIDA::IHistogram1D* getfhistoCkvRadius() {return fhistoCkvRadius; }

#endif
 
 public:

#ifdef G4ANALYSIS_USE

  AIDA::IHistogram1D* fhistoNrPhotG;
  AIDA::IHistogram1D* fhistoNBeforeMirror;  
  AIDA::IHistogram1D* fhistoWBeforeMirror;  
  AIDA::IHistogram1D* fhistoWAfterMirror;
  AIDA::IHistogram1D* fhistoCkvProdSmall;  
  AIDA::IHistogram1D* fhistoEmisZ;  
  AIDA::IHistogram1D* fhistoCkvRadius; 
  
#endif

  G4Timer* iTimer;



  G4int NumPhotBeforeAerogel;
  G4int NumPhotBeforeMirror;
  G4int NumPhotBeforeFilter;
  G4int NumPhotAfterFilter;
  G4int NumPhotAtHpd1Input;
  G4int NumPhotAtHpd2Input;
  G4int NumPhotAtHpd3Input;
  G4int NumPhotAtHpd4Input;
  G4int NumPhotAfterMirror;
  G4int NumHitTotInHpd1;
  G4int NumHitTotInHpd2;
  G4int NumHitTotInHpd3;
  G4int NumHitTotInHpd4;
  G4int NumHitInSi;
  RichTbRunConfig* rConfig;
  G4double MeanHXCoord;
  G4double MeanHYCoord;
  G4double NumHXCoord;
  G4double NumHYCoord;

};

#endif
