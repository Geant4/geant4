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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//
 // HPW to implement the choosing of an element for scattering.
#include "g4std/fstream"
#include "g4std/strstream"
#include <stdlib.h>
#include "G4HadronicProcess.hh"
//@@ add model name info, once typeinfo available #include <typeinfo.h>
 
 G4IsoParticleChange * G4HadronicProcess::theIsoResult = NULL;
 G4IsoParticleChange * G4HadronicProcess::theOldIsoResult = NULL;
 G4bool G4HadronicProcess::isoIsEnabled = true;
 void G4HadronicProcess::EnableIsotopeProductionGlobally()  {isoIsEnabled = true;}
 void G4HadronicProcess::DisableIsotopeProductionGlobally() {isoIsEnabled = false;}
 
 G4Element * G4HadronicProcess::ChooseAandZ(
  const G4DynamicParticle *aParticle, const G4Material *aMaterial )
  {
    currentZ = 0;
    currentN = 0;
    const G4int numberOfElements = aMaterial->GetNumberOfElements();
    const G4ElementVector *theElementVector = aMaterial->GetElementVector();
    
    if( numberOfElements == 1 ) 
    {
      currentZ = G4double( ((*theElementVector)[0])->GetZ());
      currentN = (*theElementVector)[0]->GetN();
      targetNucleus.SetParameters(currentN, currentZ);
      return (*theElementVector)[0];
    }
    
    const G4double *theAtomicNumberDensity = aMaterial->GetAtomicNumDensityVector();
    G4double aTemp = aMaterial->GetTemperature();
    G4double crossSectionTotal = 0;
    G4int i;
    G4std::vector<G4double> runningSum;
    for( i=0; i < numberOfElements; ++i )
    {
      runningSum.push_back(theAtomicNumberDensity[i] *
        dispatch->GetMicroscopicCrossSection( aParticle, (*theElementVector)[i], aTemp));
      crossSectionTotal+=runningSum[i];
    }
    
    G4double random = G4UniformRand();
    for( i=0; i < numberOfElements; ++i )
    { 
      if( random<=runningSum[i]/crossSectionTotal )
      {
        currentZ = G4double( ((*theElementVector)[i])->GetZ());
        currentN = ((*theElementVector)[i])->GetN();
        targetNucleus.SetParameters(currentN, currentZ);
        return (*theElementVector)[i];
      }
    }
    currentZ = G4double((*theElementVector)[numberOfElements-1]->GetZ());
    currentN = (*theElementVector)[numberOfElements-1]->GetN();
    targetNucleus.SetParameters(currentN, currentZ);
    return (*theElementVector)[numberOfElements-1];
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
    if(getenv("HadronicDoitLogging") )
    {
      G4cout << "HadronicDoitLogging "
             << GetProcessName() <<" "
             << aParticle->GetDefinition()->GetPDGEncoding()<<" "
	     << kineticEnergy<<" "
	     << aParticle->GetMomentum()<<" "
	     << targetNucleus.GetN()<<" "
	     << targetNucleus.GetZ()<<" "
	     << G4endl;
    }
    ResetNumberOfInteractionLengthLeft();
    if(isoIsOnAnyway!=-1)
    {
      if(isoIsEnabled||isoIsOnAnyway)
      {
        result = DoIsotopeCounting(result, aTrack, targetNucleus);
      }
    }
    return result;
  }

  G4VParticleChange * G4HadronicProcess::
  DoIsotopeCounting(G4VParticleChange * aResult,
                    const G4Track & aTrack,
                    const G4Nucleus & aNucleus)
  {
    // get the PC from iso-production
    if(theOldIsoResult) delete theOldIsoResult;
    if(theIsoResult) delete theIsoResult;
    theIsoResult = new G4IsoParticleChange;
    G4bool done = false;
    G4IsoResult * anIsoResult = NULL;
    for(unsigned int i=0; i<theProductionModels.size(); i++)
    {
      anIsoResult = theProductionModels[i]->GetIsotope(aTrack, aNucleus);
      if(anIsoResult!=NULL)
      {
        done = true;
        break;
      }
    }
    // if none in charge, use default iso production
    if(!done) anIsoResult = ExtractResidualNucleus(aTrack, aNucleus, aResult); 
        
    // Add all info explicitely and add typename from model called.
    theIsoResult->SetIsotope(anIsoResult->GetIsotope());
    theIsoResult->SetProductionPosition(aTrack.GetPosition());
    theIsoResult->SetProductionTime(aTrack.GetGlobalTime());
    theIsoResult->SetParentParticle(*aTrack.GetDynamicParticle());
    theIsoResult->SetMotherNucleus(anIsoResult->GetMotherNucleus());
//    theIsoResult->SetProducer(typeid(*theInteraction).name()); @@@@@@@
    G4String aWorkaround("WaitingForTypeidToBeAvailableInCompilers"); // @@@@@  workaround for DEC.
    theIsoResult->SetProducer(aWorkaround);
    
    delete anIsoResult;

    return aResult;
  }
  
  G4IsoResult * G4HadronicProcess::
  ExtractResidualNucleus(const G4Track & aTrack,
                         const G4Nucleus & aNucleus,
                         G4VParticleChange * aResult)
  {
    G4double A = aNucleus.GetN();
    G4double Z = aNucleus.GetZ();
    G4double bufferA = 0;
    G4double bufferZ = 0;
    
    // loop over aResult, and decrement A, Z accordingly
    // cash the max
    for(G4int i=0; i<aResult->GetNumberOfSecondaries(); i++)
    {
      G4Track* aSecTrack = aResult->GetSecondary(i);
      if(bufferA<aSecTrack->GetDefinition()->GetBaryonNumber())
      {
        bufferA = aSecTrack->GetDefinition()->GetBaryonNumber();
        bufferZ = aSecTrack->GetDefinition()->GetPDGCharge();
      }
      Z-=aSecTrack->GetDefinition()->GetPDGCharge();
      A-=aSecTrack->GetDefinition()->GetBaryonNumber();
    }
    
    // if the fragment was part of the final state, it is 
    // assumed to be the heaviest secondary.
    if(A<0.1)
    {
       A = bufferA;
       Z = bufferZ;
    }
    
    // prepare the IsoResult.
    char the1[100] = {""};
    G4std::ostrstream ost1(the1, 100, G4std::ios::out);
    ost1 <<Z<<"_"<<A<<"\0";
    G4String * biff = new G4String(the1);
    G4IsoResult * theResult = new G4IsoResult(*biff, aNucleus);
    
    // cleaning up.
    delete biff;
    
    return theResult;
  }


 /* end of file */
