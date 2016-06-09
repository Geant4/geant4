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
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

// 05-11-21 NeutronHP or Low Energy Prameterization Models

#include "G4NeutronHPorLElasticModel.hh"
#include "G4SystemOfUnits.hh"

G4NeutronHPorLElasticModel::G4NeutronHPorLElasticModel()
{
   //theHPElastic = new G4NeutronHPElastic();
   theHPElastic = new G4NeutronHPorLElastic();
   theLElastic = new G4LElastic();
   theHPNames = new G4NeutronHPNames();
}
G4NeutronHPorLElasticModel::~G4NeutronHPorLElasticModel()
{
   delete theHPElastic;
   delete theLElastic;
   delete theHPNames;
}

G4HadFinalState* G4NeutronHPorLElasticModel::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus)
{
   if ( aTrack.GetKineticEnergy() > 20*MeV ) 
   {
      //G4cout << "Select LE model " << G4endl;
      return theLElastic->ApplyYourself( aTrack , aTargetNucleus );
   }

   //G4int Z =  (G4int)(aTargetNucleus.GetZ()+0.5);  
   //migrate to integer A and Z
   G4int Z =  aTargetNucleus.GetZ_asInt();  
   G4String theNameOfElement = theHPNames->GetName( Z-1 ); // GetName(0) reply "Hydrogen" 

   if ( theHPElastic->IsThisElementOK( theNameOfElement ) )
   {
      //G4cout << "Select HP model " << G4endl;
      return theHPElastic->ApplyYourself( aTrack , aTargetNucleus ); 
   }
   else
   {
      //G4cout << "Select LE model " << G4endl;
      return theLElastic->ApplyYourself( aTrack , aTargetNucleus ); 
   }
}

