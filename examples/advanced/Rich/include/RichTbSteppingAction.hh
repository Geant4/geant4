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
// RichTbSteppingAction.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbSteppingAction_h
#define RichTbSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"

class RichTbAnalysisManager;
class RichTbRunConfig;
class RichTbPrimaryGeneratorAction;

class RichTbSteppingAction : public G4UserSteppingAction {

 public:

  RichTbSteppingAction(RichTbRunConfig* ,
      RichTbPrimaryGeneratorAction* );
  virtual ~RichTbSteppingAction();
  void UserSteppingAction(const G4Step* aStep);
  void RichTbGenericHisto(const G4Step* aStep);
  void  RichTbDebugHisto(const G4Step* aStep);
  RichTbAnalysisManager* getRAnalysisManager()
  {return  ranalysisManager; }
  RichTbRunConfig* getRichTbRunConfig() 
  {return richtbRunConfig; }  
  G4double getHpdPhElectronKE() 
  {return HpdPhElectronKE; }
  
private:
  RichTbAnalysisManager* ranalysisManager;
  RichTbRunConfig* richtbRunConfig;
  RichTbPrimaryGeneratorAction* rPrimGenAction;
  G4double HpdPhElectronKE;
  G4VParticleChange* uParticleChange;

};
#endif
