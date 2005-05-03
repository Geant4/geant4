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
// Rich advanced example for Geant4
// RichTbEventAction.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbEventAction_h
#define RichTbEventAction_h 1
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "RichTbRunConfig.hh"
#include "RichTbAnalysisManager.hh"

class G4Event;
class G4VVisManager;
class RichTbAnalysisManager;
class RichTbIOData;

class RichTbEventAction : public G4UserEventAction {

 public:
  RichTbEventAction();

  RichTbEventAction(RichTbRunConfig* ,
		    G4VVisManager* , RichTbIOData* );
  virtual ~RichTbEventAction();
public:
    void BeginOfEventAction(const G4Event* );
    void EndOfEventAction(const G4Event* );
    G4int GetRichCollID() {return RichTbCollID;};
  RichTbAnalysisManager* getAnalysisM()
  {return ranalysisManager; }
  G4VVisManager* getVisM()
  {return pVisManager; }
 private:
  RichTbRunConfig* runConfiguration;
  G4int RichTbCollID;
  RichTbAnalysisManager* ranalysisManager;
  G4VVisManager* pVisManager;
  RichTbIOData* rTbIOData;

};
#endif
