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
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//
//      ---------- G4NuclearStopping physics process -----
//                by Vladimir Ivanchenko, 27 April 2004
//
// Modifications:
//


// -----------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4NuclearStopping.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4ProcessManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4NuclearStopping::G4NuclearStopping(const G4String& processName, G4ProcessType aType)
  : G4VContinuousProcess(processName,aType),
    theTable("ICRU_R49"),
    model(0),
    fluctuations(true),
    initialised(false),
    factorsAreActive(false)
{
  lowEnergy        = 0.1*keV;
  highEnergy       = 10.*MeV;
  pParticleChange  = &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4NuclearStopping::~G4NuclearStopping()
{
  if(model) delete model;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NuclearStopping::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(!initialised) {
    model = new G4hNuclearStoppingModel(theTable);
    if(fluctuations) model->SetNuclearStoppingFluctuationsOn();
    else             model->SetNuclearStoppingFluctuationsOff();
    PrintInfoDefinition();
  }
  initialised = true;

  if(verboseLevel > 0) {
    G4cout << "G4NuclearStopping::BuildPhysicsTable with model <"
           << theTable << "> for " << p.GetParticleName()
	   << " fluctuations= " << fluctuations
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NuclearStopping::SetNuclearStoppingPowerModel(const G4String& dedxTable)
{
  theTable = dedxTable;
  if(initialised) {
    delete model;
    model = new G4hNuclearStoppingModel(theTable);
    if(fluctuations) model->SetNuclearStoppingFluctuationsOn();
    else             model->SetNuclearStoppingFluctuationsOff();
  }
  if(verboseLevel > 0) {
    G4cout << "G4NuclearStopping::SetNuclearStoppingPowerModel:  <"
           << theTable << ">" << G4endl;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NuclearStopping::AddSaturationFactor(const G4Material* material, G4double val)
{
  if(val > 0.0 && val < 1.0) {
    factors[material] = val;
    factorsAreActive  = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange * G4NuclearStopping::AlongStepDoIt(const G4Track& track,
                                                     const G4Step& step)
{
  fParticleChange.InitializeForAlongStep(track);
  G4double tkin = track.GetKineticEnergy();
  if(tkin > highEnergy) return &fParticleChange;

  G4double length = step.GetStepLength();
  G4double eloss  = 0.0;
  if (tkin > 0.0) eloss = length*preStepDEDX;

  if(eloss >= tkin ) {
    eloss = tkin;
    tkin = 0.0;
    if(track.GetDynamicParticle()->GetDefinition()->
             GetProcessManager()->GetAtRestProcessVector()->size())
               fParticleChange.SetStatusChange(fStopButAlive);

    else       fParticleChange.SetStatusChange(fStopAndKill);

  } else {
    tkin -= eloss;
  }
  fParticleChange.SetProposedKineticEnergy(tkin);

  if (factorsAreActive) {
     G4double factor = 1.0;
     const G4Material* material = track.GetMaterial();
     std::map<const G4Material*,G4double,std::less<const G4Material*> >::const_iterator pos;
     for (pos = factors.begin(); pos != factors.end(); pos++) {
      if((*pos).first == material) factor = (*pos).second;
     }

     eloss *= factor;
  }
  fParticleChange.SetLocalEnergyDeposit(eloss);

  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NuclearStopping::PrintInfoDefinition()
{
  if(-1 < verboseLevel) {
    G4cout << G4endl << GetProcessName() << ":   parameterization <" << theTable
           << "> is used; maxEnergy = " << highEnergy/MeV << " MeV; fluctuations= "
	   << fluctuations << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
