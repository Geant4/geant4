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
// Class Description
// Final state production model for a LEND (Low Energy Nuclear Data) 
// LEND is Geant4 interface for GIDI (General Interaction Data Interface) 
// which gives a discription of nuclear and atomic reactions, such as
//    Binary collision cross sections
//    Particle number multiplicity distributions of reaction products
//    Energy and angular distributions of reaction products
//    Derived calculational constants
// GIDI is developped at Lawrence Livermore National Laboratory
// Class Description - End

// 071025 First implementation done by T. Koi (SLAC/SCCS)
// 101118 Name modifications for release T. Koi (SLAC/PPA)

#include "G4LENDModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

double MyRNG(void*) { return  G4Random::getTheEngine()->flat(); }

G4LENDModel::G4LENDModel( G4String name )
:G4HadronicInteraction( name )
{

   proj = NULL; //will be set in an inherited class

   SetMinEnergy( 0.*eV );
   SetMaxEnergy( 20.*MeV );

   //default_evaluation = "endl99"; 
   //default_evaluation = "ENDF.B-VII.0";
   default_evaluation = "ENDF/BVII.1";

   allow_nat = false;
   allow_any = false;

   lend_manager = G4LENDManager::GetInstance();  

}

G4LENDModel::~G4LENDModel()
{
   for ( std::map< G4int , G4LENDUsedTarget* >::iterator 
         it = usedTarget_map.begin() ; it != usedTarget_map.end() ; it ++ )
   { 
      delete it->second;  
   }
}


void G4LENDModel::recreate_used_target_map()
{

   for ( std::map< G4int , G4LENDUsedTarget* >::iterator 
         it = usedTarget_map.begin() ; it != usedTarget_map.end() ; it ++ )
   { 
      delete it->second;  
   }
   usedTarget_map.clear();

   create_used_target_map();

}



void G4LENDModel::create_used_target_map()
{

   lend_manager->RequestChangeOfVerboseLevel( verboseLevel );

   size_t numberOfElements = G4Element::GetNumberOfElements();
   static const G4ElementTable* theElementTable = G4Element::GetElementTable();

   for ( size_t i = 0 ; i < numberOfElements ; ++i )
   {

      const G4Element* anElement = (*theElementTable)[i];
      G4int numberOfIsotope = anElement->GetNumberOfIsotopes(); 

      if ( numberOfIsotope > 0 )
      {
      // User Defined Abundances   
         for ( G4int i_iso = 0 ; i_iso < numberOfIsotope ; i_iso++ )
         {
            G4int iZ = anElement->GetIsotope( i_iso )->GetZ();
            G4int iA = anElement->GetIsotope( i_iso )->GetN();
            G4int iIsomer = anElement->GetIsotope( i_iso )->Getm();

            G4LENDUsedTarget* aTarget = new G4LENDUsedTarget ( proj , default_evaluation , iZ , iA , iIsomer );  
            if ( allow_nat == true ) aTarget->AllowNat();
            if ( allow_any == true ) aTarget->AllowAny();
            usedTarget_map.insert( std::pair< G4int , G4LENDUsedTarget* > ( lend_manager->GetNucleusEncoding( iZ , iA , iIsomer ) , aTarget ) );
         }
      }
      else
      {
      // Natural Abundances   
         G4NistElementBuilder* nistElementBuild = lend_manager->GetNistElementBuilder();
         G4int iZ = int ( anElement->GetZ() );
         //G4cout << nistElementBuild->GetNumberOfNistIsotopes( int ( anElement->GetZ() ) ) << G4endl;
         G4int numberOfNistIso = nistElementBuild->GetNumberOfNistIsotopes( int ( anElement->GetZ() ) ); 

         for ( G4int ii = 0 ; ii < numberOfNistIso ; ii++ )
         {
            //G4cout << nistElementBuild->GetIsotopeAbundance( iZ , nistElementBuild->GetNistFirstIsotopeN( iZ ) + i ) << G4endl;
            if ( nistElementBuild->GetIsotopeAbundance( iZ , nistElementBuild->GetNistFirstIsotopeN( iZ ) + ii ) > 0 )
            {
               G4int iMass = nistElementBuild->GetNistFirstIsotopeN( iZ ) + ii;  
               //G4cout << iZ << " " << nistElementBuild->GetNistFirstIsotopeN( iZ ) + i << " " << nistElementBuild->GetIsotopeAbundance ( iZ , iMass ) << G4endl;  
               G4int iIsomer = 0;

               G4LENDUsedTarget* aTarget = new G4LENDUsedTarget ( proj , default_evaluation , iZ , iMass );  
               if ( allow_nat == true ) aTarget->AllowNat();
               if ( allow_any == true ) aTarget->AllowAny();
               usedTarget_map.insert( std::pair< G4int , G4LENDUsedTarget* > ( lend_manager->GetNucleusEncoding( iZ , iMass , iIsomer ) , aTarget ) );

            }

         }

      }
   }



   G4cout << "Dump UsedTarget for " << GetModelName() << G4endl;
   //G4cout << "Requested Evaluation, Z , A -> Actual Evaluation, Z , A(0=Nat) , Pointer of Target" << G4endl;
   G4cout << "Requested Evaluation, Z , A -> Actual Evaluation, Z , A(0=Nat) " << G4endl;
   for ( std::map< G4int , G4LENDUsedTarget* >::iterator 
         it = usedTarget_map.begin() ; it != usedTarget_map.end() ; it ++ )
   {
      G4cout 
         << " " << it->second->GetWantedEvaluation() 
         << ", " << it->second->GetWantedZ() 
         << ", " << it->second->GetWantedA() 
         << " -> " << it->second->GetActualEvaluation() 
         << ", " << it->second->GetActualZ() 
         << ", " << it->second->GetActualA() 
         //<< ", " << it->second->GetTarget() 
         << G4endl; 
   } 

}
  

  
#include "G4IonTable.hh"
  
G4HadFinalState * G4LENDModel::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTarg )
{

   G4double temp = aTrack.GetMaterial()->GetTemperature();

   //G4int iZ = int ( aTarg.GetZ() );
   //G4int iA = int ( aTarg.GetN() );
   //migrate to integer A and Z (GetN_asInt returns number of neutrons in the nucleus since this) 
   G4int iZ = aTarg.GetZ_asInt();
   G4int iA = aTarg.GetA_asInt();
   G4int iM = 0;
   if ( aTarg.GetIsotope() != NULL ) {
      iM = aTarg.GetIsotope()->Getm();
   }

   G4double ke = aTrack.GetKineticEnergy();

   G4HadFinalState* theResult = new G4HadFinalState();

   G4GIDI_target* aTarget = usedTarget_map.find( lend_manager->GetNucleusEncoding( iZ , iA , iM ) )->second->GetTarget();

   G4double aMu = aTarget->getElasticFinalState( ke*MeV, temp, NULL, NULL );

   G4double phi = twopi*G4UniformRand();
   G4double theta = std::acos( aMu );
   //G4double sinth = std::sin( theta );

   G4ReactionProduct theNeutron( aTrack.GetDefinition() );
   theNeutron.SetMomentum( aTrack.Get4Momentum().vect() );
   theNeutron.SetKineticEnergy( ke );

   G4ParticleDefinition* pd = G4IonTable::GetIonTable()->GetIon( iZ , iA , iM );
   G4ReactionProduct theTarget( pd );

   G4double mass = pd->GetPDGMass();

// add Thermal motion 
   G4double kT = k_Boltzmann*temp;
   G4ThreeVector v ( G4RandGauss::shoot() * std::sqrt( kT*mass ) 
                   , G4RandGauss::shoot() * std::sqrt( kT*mass ) 
                   , G4RandGauss::shoot() * std::sqrt( kT*mass ) );

   theTarget.SetMomentum( v );


     G4ThreeVector the3Neutron = theNeutron.GetMomentum();
     G4double nEnergy = theNeutron.GetTotalEnergy();
     G4ThreeVector the3Target = theTarget.GetMomentum();
     G4double tEnergy = theTarget.GetTotalEnergy();
     G4ReactionProduct theCMS;
     G4double totE = nEnergy+tEnergy;
     G4ThreeVector the3CMS = the3Target+the3Neutron;
     theCMS.SetMomentum(the3CMS);
     G4double cmsMom = std::sqrt(the3CMS*the3CMS);
     G4double sqrts = std::sqrt((totE-cmsMom)*(totE+cmsMom));
     theCMS.SetMass(sqrts);
     theCMS.SetTotalEnergy(totE);

       theNeutron.Lorentz(theNeutron, theCMS);
       theTarget.Lorentz(theTarget, theCMS);
       G4double en = theNeutron.GetTotalMomentum(); // already in CMS.
       G4ThreeVector cms3Mom=theNeutron.GetMomentum(); // for neutron direction in CMS
       G4double cms_theta=cms3Mom.theta();
       G4double cms_phi=cms3Mom.phi();
       G4ThreeVector tempVector;
       tempVector.setX(std::cos(theta)*std::sin(cms_theta)*std::cos(cms_phi)
                       +std::sin(theta)*std::cos(phi)*std::cos(cms_theta)*std::cos(cms_phi)
                       -std::sin(theta)*std::sin(phi)*std::sin(cms_phi)  );
       tempVector.setY(std::cos(theta)*std::sin(cms_theta)*std::sin(cms_phi)
                       +std::sin(theta)*std::cos(phi)*std::cos(cms_theta)*std::sin(cms_phi)
                       +std::sin(theta)*std::sin(phi)*std::cos(cms_phi)  );
       tempVector.setZ(std::cos(theta)*std::cos(cms_theta)
                       -std::sin(theta)*std::cos(phi)*std::sin(cms_theta)  );
       tempVector *= en;
       theNeutron.SetMomentum(tempVector);
       theTarget.SetMomentum(-tempVector);
       G4double tP = theTarget.GetTotalMomentum();
       G4double tM = theTarget.GetMass();
       theTarget.SetTotalEnergy(std::sqrt((tP+tM)*(tP+tM)-2.*tP*tM));
       theNeutron.Lorentz(theNeutron, -1.*theCMS);
       theTarget.Lorentz(theTarget, -1.*theCMS);

     theResult->SetEnergyChange(theNeutron.GetKineticEnergy());
     theResult->SetMomentumChange(theNeutron.GetMomentum().unit());
     G4DynamicParticle* theRecoil = new G4DynamicParticle;

     theRecoil->SetDefinition( G4IonTable::GetIonTable()->GetIon( iZ , iA , iM , iZ ) );
     theRecoil->SetMomentum( theTarget.GetMomentum() );

     theResult->AddSecondary( theRecoil );

   return theResult; 

}
