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
// $Id: Tst50PositronPenelope.cc,v 1.1 2004-06-02 10:16:09 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria.Grazia.Pia@cern.ch
//
// History:
// -----------
// 22 Feb 2003 MGP          Designed for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst50PositronPenelope.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4PenelopeIonisation.hh"
#include "G4PenelopeBremsstrahlung.hh"
#include "G4PenelopeAnnihilation.hh"

Tst50PositronPenelope::Tst50PositronPenelope(const G4String& name): 
G4VPhysicsConstructor(name)
{ }

Tst50PositronPenelope::~Tst50PositronPenelope()
{ }

void Tst50PositronPenelope::ConstructProcess()
{
  
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "e+") 
	{
	  manager -> AddProcess(new G4MultipleScattering,    -1, 1,1);
          manager -> AddProcess(new G4PenelopeIonisation,    -1, 2,2);
	  manager -> AddProcess(new G4PenelopeBremsstrahlung,-1,-1,3);  
          manager -> AddProcess(new G4PenelopeAnnihilation,   0,-1,4);
	}   
    }
}
