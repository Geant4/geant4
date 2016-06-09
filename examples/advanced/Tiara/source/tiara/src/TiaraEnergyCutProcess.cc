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
// $Id: TiaraEnergyCutProcess.cc,v 1.4 2003/11/07 12:39:11 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00 $
//

#include "TiaraEnergyCutProcess.hh"

TiaraEnergyCutProcess::TiaraEnergyCutProcess(G4double minEnergyCut)
 : 
  G4VProcess("TiaraEnergyCutProcess"), 
  fMinEnergyCut(minEnergyCut)
{
  G4VProcess::pParticleChange = new G4ParticleChange;
  if (!G4VProcess::pParticleChange) {
    G4Exception("ERROR:TiaraEnergyCutProcess::TiaraEnergyCutProcess new failed to create G4ParticleChange!");
  }
}

TiaraEnergyCutProcess::~TiaraEnergyCutProcess()
{
  delete pParticleChange;
}

G4double TiaraEnergyCutProcess::
PostStepGetPhysicalInteractionLength(const G4Track&,
				     G4double,
				     G4ForceCondition* condition)
{
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange * 
TiaraEnergyCutProcess::PostStepDoIt(const G4Track& aTrack, const G4Step &)
{
  pParticleChange->Initialize(aTrack);

  if (aTrack.GetKineticEnergy() < fMinEnergyCut) {
    pParticleChange->SetStatusChange(fStopAndKill);
  }
  
  return G4VProcess::pParticleChange;
}

const G4String &TiaraEnergyCutProcess::GetName() const {
  return theProcessName;
}


G4double TiaraEnergyCutProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
				      G4double  ,
				      G4double  ,
				      G4double& ,
				      G4GPILSelection*) {
  return -1.0;
}

G4double TiaraEnergyCutProcess::
AtRestGetPhysicalInteractionLength(const G4Track&,
				   G4ForceCondition*) {
  return -1.0;
}

G4VParticleChange* TiaraEnergyCutProcess::AtRestDoIt(const G4Track&,
					       const G4Step&) {
  return 0;
}

G4VParticleChange* TiaraEnergyCutProcess::AlongStepDoIt(const G4Track&,
						  const G4Step&) {
  return 0;
}

