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
//
// $Id: G4EnergyRangeManager.cc,v 1.16 2010-11-22 07:45:45 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Energy Range Manager
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 22-Nov-1996
 // Last modified: 24-Mar-1997
 // fix in the counter-hndling: H.P. Wellisch 04-Apr-97
 // throw an exception if no model found:  J.L. Chuma  04-Apr-97
 
#include "G4EnergyRangeManager.hh"
#include "Randomize.hh"
#include "G4HadronicException.hh"


G4EnergyRangeManager::G4EnergyRangeManager()
 : theHadronicInteractionCounter(0)
{
  for (G4int i = 0; i < G4EnergyRangeManager::MAX_NUMBER_OF_MODELS; i++)
    theHadronicInteraction[i] = 0;
}


G4EnergyRangeManager::G4EnergyRangeManager(const G4EnergyRangeManager& right)
{
  if (this != &right) {
    for (G4int i=0; i<theHadronicInteractionCounter; ++i)
      theHadronicInteraction[i] = right.theHadronicInteraction[i];
    theHadronicInteractionCounter = right.theHadronicInteractionCounter;
  }
}
 

G4EnergyRangeManager&G4EnergyRangeManager::operator=(
   const G4EnergyRangeManager &right )
{
  if (this != &right) {
    for (G4int i=0; i<theHadronicInteractionCounter; ++i)
      theHadronicInteraction[i] = right.theHadronicInteraction[i];
    theHadronicInteractionCounter = right.theHadronicInteractionCounter;
  }
  return *this;
}


void G4EnergyRangeManager::RegisterMe(G4HadronicInteraction* a)
{
    if( theHadronicInteractionCounter+1 > MAX_NUMBER_OF_MODELS )
    {
      throw G4HadronicException(__FILE__, __LINE__,"RegisterMe: TOO MANY MODELS");
    }
    theHadronicInteraction[ theHadronicInteractionCounter++ ] = a;
}

 
G4HadronicInteraction*
  G4EnergyRangeManager::GetHadronicInteraction(
   const G4double kineticEnergy,
   const G4Material *aMaterial,
   const G4Element *anElement ) const
  {
    G4int counter = GetHadronicInteractionCounter();
    if( counter == 0 )
      throw G4HadronicException(__FILE__, __LINE__,
                               "GetHadronicInteraction: NO MODELS STORED");

    G4int cou = 0, memory = 0, memor2 = 0;
    G4double emi1 = 0.0, ema1 = 0.0, emi2 = 0.0, ema2 = 0.0;
    for( G4int i=0; i<counter; i++ )
    {
      G4double low  = theHadronicInteraction[i]->GetMinEnergy( aMaterial, anElement );
      // Work-around for particles with 0 kinetic energy, which still
      // require a model to return a ParticleChange
      if (low == 0.) low = -DBL_MIN;
      G4double high = theHadronicInteraction[i]->GetMaxEnergy( aMaterial, anElement );
      if( low < kineticEnergy && high >= kineticEnergy )
      {
        ++cou;
        emi2 = emi1;
        ema2 = ema1;
        emi1 = low;
        ema1 = high;
        memor2 = memory;
        memory = i;
      }
    }
    G4int m=-1;
    G4double rand;
    switch ( cou )
    {
     case 0:
       G4cout<<"G4EnergyRangeManager:GetHadronicInteraction: counter="<<counter<<", Ek="
             <<kineticEnergy<<", Material = "<<aMaterial->GetName()<<", Element = "
             <<anElement->GetName()<<G4endl;
       for( G4int j=0; j<counter; j++ )
       {
         G4HadronicInteraction* HInt=theHadronicInteraction[j];
         G4cout<<"*"<<j<<"* low=" <<HInt->GetMinEnergy(aMaterial,anElement)
               <<", high="<<HInt->GetMaxEnergy(aMaterial,anElement)<<G4endl;
       }
       throw G4HadronicException(__FILE__, __LINE__,
          "GetHadronicInteraction: No Model found");
       return 0;
     case 1:
       m = memory;
       break;
     case 2:
       if( (emi2<=emi1 && ema2>=ema1) || (emi2>=emi1 && ema2<=ema1) )
       {
         G4cout<<"G4EnergyRangeManager:GetHadronicInteraction: counter="<<counter<<", Ek="
               <<kineticEnergy<<", Material = "<<aMaterial->GetName()<<", Element = "
               <<anElement->GetName()<<G4endl;
         if(counter) for( G4int j=0; j<counter; j++ )
         {
           G4HadronicInteraction* HInt=theHadronicInteraction[j];
           G4cout<<"*"<<j<<"* low=" <<HInt->GetMinEnergy(aMaterial,anElement)
               <<", high="<<HInt->GetMaxEnergy(aMaterial,anElement)<<G4endl;
         }
         throw G4HadronicException(__FILE__, __LINE__,
               "GetHadronicInteraction: Energy ranges of two models fully overlapping");
       }
       rand = G4UniformRand();
       if( emi1 < emi2 )
       {
         if( (ema1-kineticEnergy)/(ema1-emi2)<rand )
           m = memor2;
         else
           m = memory;
       } else {
         if( (ema2-kineticEnergy)/(ema2-emi1)<rand )
           m = memory;
         else
           m = memor2;
       }
       break;
     default:
      throw G4HadronicException(__FILE__, __LINE__,
        "GetHadronicInteraction: More than two competing models in this energy range");
    }
    return theHadronicInteraction[m];
  } 
 
 /* end of file */
 
