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
// $Id: G4MassImportanceProcess.hh,v 1.2 2002-04-09 17:40:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MassImportanceProcess
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassImportanceProcess_hh
#define G4MassImportanceProcess_hh

#include "G4VProcess.hh"
#include "G4ImportancePostStepDoIt.hh"

class G4VImportanceAlgorithm;
class G4ImportanceFinder;
class G4VIStore;

class G4MassImportanceProcess : public G4VProcess
{

public:  // with description

  G4MassImportanceProcess(const G4VImportanceAlgorithm &aImportanceAlgorithm,
			  const G4VIStore &aIstore,
			  const G4String &aName = "MassImportanceProcess");
  ~G4MassImportanceProcess();
  
  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double   previousStepSize,
				       G4ForceCondition* condition);
  
  virtual G4VParticleChange * PostStepDoIt(const G4Track&, const G4Step&);

public:  // without description

  //  no operation in  AtRestDoIt and  AlongStepDoIt

  G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
					G4double  ,
					G4double  ,
					G4double& ,
					G4GPILSelection*) {return -1.0;}
  
  G4double AtRestGetPhysicalInteractionLength(const G4Track& ,
				     G4ForceCondition*) {return -1.0;}
  
  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) {return 0;}
  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) {return 0;}
  
private:

  G4MassImportanceProcess(const G4MassImportanceProcess &);
  G4MassImportanceProcess &operator=(const G4MassImportanceProcess &);

private:

  G4ParticleChange *fParticleChange;
  G4ImportancePostStepDoIt fImportancePostStepDoIt;
  const G4VImportanceAlgorithm &fImportanceAlgorithm;
  G4ImportanceFinder *fImportanceFinder;
};

#endif
