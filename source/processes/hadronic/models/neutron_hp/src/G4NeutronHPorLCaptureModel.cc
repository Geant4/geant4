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

