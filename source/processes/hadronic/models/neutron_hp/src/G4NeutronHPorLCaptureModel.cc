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
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

#include "G4NeutronHPorLCaptureModel.hh"

G4NeutronHPorLCaptureModel::G4NeutronHPorLCaptureModel()
{
   theHPCapture = new G4NeutronHPorLCapture();
   theLCapture = new G4LCapture();
   theHPNames = new G4NeutronHPNames();
}
G4NeutronHPorLCaptureModel::~G4NeutronHPorLCaptureModel()
{
   delete theHPCapture;
   delete theLCapture;
   delete theHPNames;
}

G4HadFinalState* G4NeutronHPorLCaptureModel::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus)
{
   if ( aTrack.GetKineticEnergy() > 20*MeV ) 
   {
      //G4cout << "Select LE model " << G4endl;
      return theLCapture->ApplyYourself( aTrack , aTargetNucleus );
   }

   G4int Z =  (G4int)(aTargetNucleus.GetZ()+0.5);  
   G4String theNameOfElement = theHPNames->GetName( Z-1 ); // GetName(0) reply "Hydrogen" 

   if ( theHPCapture->IsThisElementOK( theNameOfElement ) )
   {
      //G4cout << "Select HP model " << G4endl;
      return theHPCapture->ApplyYourself( aTrack , aTargetNucleus ); 
   }
   else
   {
      //G4cout << "Select LE model " << G4endl;
      return theLCapture->ApplyYourself( aTrack , aTargetNucleus ); 
   }
}

