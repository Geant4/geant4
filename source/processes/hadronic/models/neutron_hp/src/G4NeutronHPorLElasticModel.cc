//
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

// 05-11-21 NeutronHP or Low Energy Prameterization Models

#include "G4NeutronHPorLElasticModel.hh"

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

   G4int Z =  (G4int)(aTargetNucleus.GetZ()+0.5);  
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

