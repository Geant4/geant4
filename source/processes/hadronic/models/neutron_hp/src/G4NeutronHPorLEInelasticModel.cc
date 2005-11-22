//
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

#include "G4NeutronHPorLEInelasticModel.hh"

G4NeutronHPorLEInelasticModel::G4NeutronHPorLEInelasticModel()
{
   //theHPElastic = new G4NeutronHPElastic();
   theHPInelastic = new G4NeutronHPorLEInelastic();
   theLEInelastic = new G4LENeutronInelastic();
   theHPNames = new G4NeutronHPNames();
}
G4NeutronHPorLEInelasticModel::~G4NeutronHPorLEInelasticModel()
{
   delete theHPInelastic;
   delete theLEInelastic;
   delete theHPNames;
}

G4HadFinalState* G4NeutronHPorLEInelasticModel::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus)
{
   if ( aTrack.GetKineticEnergy() > 20*MeV ) 
   {
      //G4cout << "Select LE model " << G4endl;
      return theLEInelastic->ApplyYourself( aTrack , aTargetNucleus );
   }

   G4int Z =  (G4int)(aTargetNucleus.GetZ()+0.5);  
   G4String theNameOfElement = theHPNames->GetName( Z-1 ); // GetName(0) reply "Hydrogen" 

   if ( theHPInelastic->IsThisElementOK( theNameOfElement ) )
   {
      //G4cout << "Select HP model " << G4endl;
      return theHPInelastic->ApplyYourself( aTrack , aTargetNucleus ); 
   }
   else
   {
      //G4cout << "Select LE model " << G4endl;
      return theLEInelastic->ApplyYourself( aTrack , aTargetNucleus ); 
   }
}

