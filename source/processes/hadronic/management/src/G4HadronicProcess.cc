// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronicProcess.cc,v 1.1 1999-01-07 16:11:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // HPW to implement the choosing of an element for scattering.
#include "G4HadronicProcess.hh"
 
 G4Element * G4HadronicProcess::ChooseAandZ(
  const G4DynamicParticle *aParticle, const G4Material *aMaterial )
  {
    currentZ = 0;
    currentN = 0;
    const G4int numberOfElements = aMaterial->GetNumberOfElements();
    const G4ElementVector *theElementVector = aMaterial->GetElementVector();
    
    if( numberOfElements == 1 ) 
    {
      currentZ = G4double((*theElementVector)(0)->GetZ());
      currentN = (*theElementVector)(0)->GetN();
      targetNucleus.SetParameters(currentN, currentZ);
      return (*theElementVector)(0);
    }
    
    const G4double *theAtomicNumberDensity = aMaterial->GetAtomicNumDensityVector();
    G4double crossSectionTotal = 0;
    G4int i;
    for( i=0; i < numberOfElements; ++i )
      crossSectionTotal += theAtomicNumberDensity[i] *
        dispatch->GetMicroscopicCrossSection( aParticle, (*theElementVector)(i) );
    
    G4double crossSectionSum= 0.;
    G4double random = G4UniformRand()*crossSectionTotal;
    for( i=0; i < numberOfElements; ++i )
    { 
      crossSectionSum += theAtomicNumberDensity[i] *
        dispatch->GetMicroscopicCrossSection( aParticle, (*theElementVector)(i) );
      if( random<=crossSectionSum )
      {
        currentZ = G4double((*theElementVector)(i)->GetZ());
        currentN = (*theElementVector)(i)->GetN();
        targetNucleus.SetParameters(currentN, currentZ);
        return (*theElementVector)(i);
      }
    }
    currentZ = G4double((*theElementVector)(numberOfElements-1)->GetZ());
    currentN = (*theElementVector)(numberOfElements-1)->GetN();
    targetNucleus.SetParameters(currentN, currentZ);
    return (*theElementVector)(numberOfElements-1);
  }
 
 G4VParticleChange *G4HadronicProcess::GeneralPostStepDoIt(
  const G4Track &aTrack, const G4Step &aStep )
  {
    const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
    G4Material *aMaterial = aTrack.GetMaterial();
    G4double kineticEnergy = aParticle->GetKineticEnergy();
    G4Element * anElement = ChooseAandZ( aParticle, aMaterial );
    theInteraction = ChooseHadronicInteraction( kineticEnergy,
                                                          aMaterial, anElement );
    G4VParticleChange *result =
      theInteraction->ApplyYourself( aTrack, targetNucleus);
    ResetNumberOfInteractionLengthLeft();
    return result;
  }

 /* end of file */
