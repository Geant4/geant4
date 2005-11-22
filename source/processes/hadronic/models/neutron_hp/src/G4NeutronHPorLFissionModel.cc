//
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

#include "G4NeutronHPorLFissionModel.hh"

G4NeutronHPorLFissionModel::G4NeutronHPorLFissionModel()
{
   theHPFission = new G4NeutronHPorLFission();
   theLFission = new G4LFission();
   theHPNames = new G4NeutronHPNames();
}
G4NeutronHPorLFissionModel::~G4NeutronHPorLFissionModel()
{
   delete theHPFission;
   delete theLFission;
   delete theHPNames;
}

G4HadFinalState* G4NeutronHPorLFissionModel::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus)
{
   if ( aTrack.GetKineticEnergy() > 20*MeV ) 
   {
      //G4cout << "Select LE model " << G4endl;
      return theLFission->ApplyYourself( aTrack , aTargetNucleus );
   }

   G4int Z =  (G4int)(aTargetNucleus.GetZ()+0.5);  
   G4String theNameOfElement = theHPNames->GetName( Z-1 ); // GetName(0) reply "Hydrogen" 

   if ( theHPFission->IsThisElementOK( theNameOfElement ) )
   {
      //G4cout << "Select HP model " << G4endl;
      return theHPFission->ApplyYourself( aTrack , aTargetNucleus ); 
   }
   else
   {
      //G4cout << "Select LE model " << G4endl;
      return theLFission->ApplyYourself( aTrack , aTargetNucleus ); 
   }
}

