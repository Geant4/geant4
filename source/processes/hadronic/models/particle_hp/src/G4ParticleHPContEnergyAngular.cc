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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 080721 To be "ClearHistories" effective, the selection scheme of angular distribution is modified by T. Koi
//
// P. Arce, Dec-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPContEnergyAngular.hh"

G4ReactionProduct*
G4ParticleHPContEnergyAngular::Sample(G4double anEnergy, G4double massCode,
                                      G4double /*mass*/)
{
  G4ReactionProduct* result;
  G4int i(0);
  G4int it(0);
  for (i = 0; i < nEnergy; i++) {
    it = i;
#ifdef PHP_AS_HP 
    if(theAngular[i].GetEnergy() > anEnergy) break;
#else
    if(theAngular[i].GetEnergy() >= anEnergy) break;
#endif
  }

  G4double targetMass = GetTarget()->GetMass();
  if (it == 0) {
    theAngular[0].SetTarget(GetTarget());
    theAngular[0].SetTargetCode(theTargetCode);
    theAngular[0].SetPrimary(GetProjectileRP());
    result = theAngular[0].Sample(anEnergy, massCode, targetMass, 
                                  theAngularRep, theInterpolation);
    currentMeanEnergy.Put(theAngular[0].MeanEnergyOfThisInteraction() );

  } else {
     // interpolation through alternating sampling. This needs improvement @@@
     // This is the cause of the He3 problem !!!!!!!!
     // See to it, if you can improve this.
     //G4double random = G4UniformRand();
     //G4double deltaE = theAngular[it].GetEnergy()-theAngular[it-1].GetEnergy();
     //G4double offset = theAngular[it].GetEnergy()-anEnergy;
     //if(random<offset/deltaE) it--;
     //--- create new 
     //  if (theManager.GetScheme(0) != LINLIN) {
     //    // asserted in G4ParticleHPContEnergyAngular::init there is only one range
#ifdef PHP_AS_HP
    theAngular[it].SetTarget(GetTarget());
    theAngular[it].SetTargetCode(theTargetCode);
    theAngular[it].SetPrimary(GetProjectileRP());
    result = theAngular[it].Sample(anEnergy, massCode, targetMass, 
				    theAngularRep, theInterpolation);
    currentMeanEnergy.Put(theAngular[it].MeanEnergyOfThisInteraction() ); 
#else
    if (fCacheAngular.Get() == NULL) {
      G4ParticleHPContAngularPar* angpar = new G4ParticleHPContAngularPar(theProjectile);
      fCacheAngular.Put(angpar);
    }
    fCacheAngular.Get()->SetInterpolation(theInterpolation);
    fCacheAngular.Get()->BuildByInterpolation(anEnergy, theManager.GetScheme(0),
                                              (theAngular[it-1]), (theAngular[it]) );
    fCacheAngular.Get()->SetTarget(GetTarget());
    fCacheAngular.Get()->SetTargetCode(theTargetCode);
    fCacheAngular.Get()->SetPrimary(GetProjectileRP());

    result = fCacheAngular.Get()->Sample(anEnergy, massCode, targetMass, 
                                         theAngularRep, theInterpolation);
    currentMeanEnergy.Put(fCacheAngular.Get()->MeanEnergyOfThisInteraction() );
#endif
  } // end (it != 0) branch  

  return result;
}


G4double G4ParticleHPContEnergyAngular::MeanEnergyOfThisInteraction()
{
  G4double result(0);
  if (currentMeanEnergy.Get() < -1) {
     throw G4HadronicException(__FILE__, __LINE__,
                               "G4ParticleHPContEnergyAngular: Logical error in Product class");
  } else {
    result = currentMeanEnergy.Get();
  }
  currentMeanEnergy.Put(-2);
  return result;
}


void G4ParticleHPContEnergyAngular::ClearHistories()
{ 
  if (theAngular!= NULL) { 
    for (G4int i = 0; i < nEnergy; i++) theAngular[i].ClearHistories(); 
  }

  // Added fCacheAngular ClearHistories() - this is the one actually used!
  // Maybe theAngular does not even need ClearHistories()?
  if (fCacheAngular.Get() != 0) fCacheAngular.Get()->ClearHistories();
}

