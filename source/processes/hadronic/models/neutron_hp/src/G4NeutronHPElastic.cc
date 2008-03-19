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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 070523 bug fix for G4FPE_DEBUG on by A. Howard ( and T. Koi)
// 080319 Compilation warnings - gcc-4.3.0 fix by T. Koi
//
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticFS.hh"

  G4NeutronHPElastic::G4NeutronHPElastic()
    :G4HadronicInteraction("NeutronHPElastic")
  {
    overrideSuspension = false;
    G4NeutronHPElasticFS * theFS = new G4NeutronHPElasticFS;
    if(!getenv("G4NEUTRONHPDATA")) 
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
    dirName = getenv("G4NEUTRONHPDATA");
    G4String tString = "/Elastic/";
    dirName = dirName + tString;
//    G4cout <<"G4NeutronHPElastic::G4NeutronHPElastic testit "<<dirName<<G4endl;
    numEle = G4Element::GetNumberOfElements();
    theElastic = new G4NeutronHPChannel[numEle];
    for (G4int i=0; i<numEle; i++)
    {
      theElastic[i].Init((*(G4Element::GetElementTable()))[i], dirName);
      while(!theElastic[i].Register(theFS)) ;
    }
    delete theFS;
    SetMinEnergy(0.*eV);
    SetMaxEnergy(20.*MeV);
  }
  
  G4NeutronHPElastic::~G4NeutronHPElastic()
  {
    delete [] theElastic;
  }
  
  #include "G4NeutronHPThermalBoost.hh"
  
  G4HadFinalState * G4NeutronHPElastic::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& )
  {
    const G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = theMaterial->GetNumberOfElements();
    G4int index = theMaterial->GetElement(0)->GetIndex();
    if(n!=1)
    {
      G4int i;
      xSec = new G4double[n];
      G4double sum=0;
      const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
      G4double rWeight;    
      G4NeutronHPThermalBoost aThermalE;
      for (i=0; i<n; i++)
      {
        index = theMaterial->GetElement(i)->GetIndex();
        rWeight = NumAtomsPerVolume[i];
        xSec[i] = theElastic[index].GetXsec(aThermalE.GetThermalEnergy(aTrack,
  		                                                     theMaterial->GetElement(i),
  								     theMaterial->GetTemperature()));
        xSec[i] *= rWeight;
        sum+=xSec[i];
      }
      G4double random = G4UniformRand();
      G4double running = 0;
      for (i=0; i<n; i++)
      {
        running += xSec[i];
        index = theMaterial->GetElement(i)->GetIndex();
        //if(random<=running/sum) break;
        if( sum == 0 || random <= running/sum ) break;
      }
      delete [] xSec;
      // it is element-wise initialised.
    }
    G4HadFinalState* finalState = theElastic[index].ApplyYourself(aTrack);
    if (overrideSuspension) finalState->SetStatusChange(isAlive);
    return finalState; 
  }
