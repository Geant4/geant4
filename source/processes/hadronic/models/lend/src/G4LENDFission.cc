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
#include "G4LENDFission.hh"
#include "G4SystemOfUnits.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"
  
G4HadFinalState * G4LENDFission::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg )
{

   G4double temp = aTrack.GetMaterial()->GetTemperature();

   //migrate to integer A and Z (GetN_asInt returns number of neutrons in the nucleus since this) 
   G4int iZ = aTarg.GetZ_asInt();
   G4int iA = aTarg.GetA_asInt();
   //G4int iM = aTarg.GetM_asInt();
   G4int iM = 0;
   if ( aTarg.GetIsotope() != NULL ) {
      iM = aTarg.GetIsotope()->Getm();
   }

   G4double ke = aTrack.GetKineticEnergy();

   G4HadFinalState* theResult = &theParticleChange;
   theResult->Clear();

   G4GIDI_target* aTarget = get_target_from_map( lend_manager->GetNucleusEncoding( iZ , iA , iM ) );
   if ( aTarget == NULL ) return returnUnchanged( aTrack , theResult );
   std::vector<G4GIDI_Product>* products = aTarget->getFissionFinalState( ke*MeV, temp, MyRNG, NULL );
   if ( products != NULL ) 
   {
      for ( G4int j = 0; j < int( products->size() ); j++ ) 
      {
         G4int jZ = (*products)[j].Z; 
         G4int jA = (*products)[j].A; 
         G4int jM = (*products)[j].m; 

         //G4cout << "Z = "    << (*products)[j].Z 
         //       << ", A = "  << (*products)[j].A 
         //       << ", EK = " << (*products)[j].kineticEnergy << " [MeV]" 
         //       << ", px = " << (*products)[j].px
         //       << ", py = " << (*products)[j].py
         //       << ", pz = " << (*products)[j].pz
         //       << ", birthTimeSec = " << (*products)[j].birthTimeSec << " [second]" 
         //       << G4endl;

         G4DynamicParticle* theSec = new G4DynamicParticle;

         if ( jZ > 0 )
         {
            theSec->SetDefinition( G4IonTable::GetIonTable()->GetIon( jZ, jA , jM ) );
         } 
         else if ( jA == 1 && jZ == 0 )
         {
            theSec->SetDefinition( G4Neutron::Neutron() );
         } 
         else
         {
            theSec->SetDefinition( G4Gamma::Gamma() );
         } 

         theSec->SetMomentum( G4ThreeVector( (*products)[j].px*MeV , (*products)[j].py*MeV , (*products)[j].pz*MeV ) );
         //G4cout << theSec->GetDefinition()->GetParticleName() << G4endl;
         theResult->AddSecondary( theSec );
         //Set time for delayed neutrons
         //Current implementation is a little tricky, 
         if ( (*products)[j].birthTimeSec != 0 ) {
            G4double time = (*products)[j].birthTimeSec*second + aTrack.GetGlobalTime();
            theResult->GetSecondary(theResult->GetNumberOfSecondaries()-1)->SetTime(time);
         }
      } 
   }
   delete products;

   theResult->SetStatusChange( stopAndKill );

   return theResult; 

}
const std::pair<G4double, G4double> G4LENDFission::GetFatalEnergyCheckLevels() const
{
        // max energy non-conservation is mass of heavy nucleus
        //return std::pair<G4double, G4double>(5*perCent,250*GeV);
        return std::pair<G4double, G4double>(5*perCent,DBL_MAX);
}
