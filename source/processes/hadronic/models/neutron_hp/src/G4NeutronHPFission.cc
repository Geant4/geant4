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
// 08-08-06 delete unnecessary and harmed declaration; Bug Report[857]
//
#include "G4NeutronHPFission.hh"
#include "G4SystemOfUnits.hh"

#include "G4NeutronHPManager.hh"

  G4NeutronHPFission::G4NeutronHPFission()
    :G4HadronicInteraction("NeutronHPFission")
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 20.*MeV );
    if(!getenv("G4NEUTRONHPDATA")) 
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
    dirName = getenv("G4NEUTRONHPDATA");
    G4String tString = "/Fission";
    dirName = dirName + tString;
    numEle = G4Element::GetNumberOfElements();
    //theFission = new G4NeutronHPChannel[numEle];

    //for (G4int i=0; i<numEle; i++)
    //{ 
      //if((*(G4Element::GetElementTable()))[i]->GetZ()>89)
    //  if((*(G4Element::GetElementTable()))[i]->GetZ()>87) //TK modified for ENDF-VII
    //  {
    //    theFission[i].Init((*(G4Element::GetElementTable()))[i], dirName);
    //    theFission[i].Register(&theFS);
    //  }
    //}

    for ( G4int i = 0 ; i < numEle ; i++ ) 
    {
      theFission.push_back( new G4NeutronHPChannel );
      if((*(G4Element::GetElementTable()))[i]->GetZ()>87) //TK modified for ENDF-VII
      {
       (*theFission[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
       (*theFission[i]).Register(&theFS);
      }
    }
  }
  
  G4NeutronHPFission::~G4NeutronHPFission()
  {
    //delete [] theFission;
     for ( std::vector<G4NeutronHPChannel*>::iterator 
           it = theFission.begin() ; it != theFission.end() ; it++ )
     {
        delete *it;
     }
     theFission.clear();
  }
  
  #include "G4NeutronHPThermalBoost.hh"
  G4HadFinalState * G4NeutronHPFission::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus )
  {

    if ( numEle < (G4int)G4Element::GetNumberOfElements() ) addChannelForNewElement();

    G4NeutronHPManager::GetInstance()->OpenReactionWhiteBoard();
    const G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = theMaterial->GetNumberOfElements();
    G4int index = theMaterial->GetElement(0)->GetIndex();
    if(n!=1)
    {
      xSec = new G4double[n];
      G4double sum=0;
      G4int i;
      const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
      G4double rWeight;    
      G4NeutronHPThermalBoost aThermalE;
      for (i=0; i<n; i++)
      {
        index = theMaterial->GetElement(i)->GetIndex();
        rWeight = NumAtomsPerVolume[i];
        xSec[i] = (*theFission[index]).GetXsec(aThermalE.GetThermalEnergy(aTrack,
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
        if( sum == 0 ||  random <= running/sum ) break;
      }
      delete [] xSec;
    }
    //return theFission[index].ApplyYourself(aTrack);
    G4HadFinalState* result = (*theFission[index]).ApplyYourself(aTrack);
    aNucleus.SetParameters(G4NeutronHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA(),G4NeutronHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargZ());
    G4NeutronHPManager::GetInstance()->CloseReactionWhiteBoard();
    return result; 
  }

const std::pair<G4double, G4double> G4NeutronHPFission::GetFatalEnergyCheckLevels() const
{
        // max energy non-conservation is mass of heavy nucleus
        //return std::pair<G4double, G4double>(5*perCent,250*GeV);
        return std::pair<G4double, G4double>(5*perCent,DBL_MAX);
}



void G4NeutronHPFission::addChannelForNewElement()
{
   for ( G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; i++ ) 
   {
      theFission.push_back( new G4NeutronHPChannel );
      if ( (*(G4Element::GetElementTable()))[i]->GetZ() > 87 ) //TK modified for ENDF-VII
      {
         G4cout << "G4NeutronHPFission Prepairing Data for the new element of " << (*(G4Element::GetElementTable()))[i]->GetName() << G4endl;
         (*theFission[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
         (*theFission[i]).Register(&theFS);
      }
   }
   numEle = (G4int)G4Element::GetNumberOfElements();
}
