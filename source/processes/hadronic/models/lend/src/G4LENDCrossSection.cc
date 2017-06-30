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
// Cross Section for LEND (Low Energy Nuclear Data)
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

#include "G4LENDCrossSection.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4ElementTable.hh"
#include "G4HadronicException.hh"

//TK110811
G4bool G4LENDCrossSection::IsIsoApplicable( const G4DynamicParticle* dp, G4int /*iZ*/ , G4int /*aA*/ , 
                                            const G4Element* /*element*/ , const G4Material* /*material*/ )
{
   G4double eKin = dp->GetKineticEnergy();
   if ( dp->GetDefinition() != proj ) return false;
   if ( eKin > GetMaxKinEnergy() || eKin < GetMinKinEnergy() ) return false;

   return true;
}

G4double G4LENDCrossSection::GetIsoCrossSection( const G4DynamicParticle* dp , G4int iZ , G4int iA ,
                                                 const G4Isotope* isotope , const G4Element* /*elment*/ , const G4Material* material )
{

   G4double xs = 0.0;
   G4double ke = dp->GetKineticEnergy();
   G4double temp = material->GetTemperature();
   G4int iM = isotope->Getm();

   G4GIDI_target* aTarget = usedTarget_map.find( lend_manager->GetNucleusEncoding( iZ , iA , iM ) )->second->GetTarget();

   xs = getLENDCrossSection ( aTarget , ke , temp );

   return xs;
}


/*
G4bool G4LENDCrossSection::IsApplicable(const G4DynamicParticle*aP, const G4Element*)
{
   G4bool result = true;
   G4double eKin = aP->GetKineticEnergy();
   if( eKin > GetMaxKinEnergy() || aP->GetDefinition() != proj ) result = false;
   return result;
}
*/

G4LENDCrossSection::G4LENDCrossSection( const G4String nam )
:G4VCrossSectionDataSet( nam )
{

   proj = NULL; //will be set in an inherited class
   //default_evaluation = "endl99";
   //default_evaluation = "ENDF.B-VII.0";
   default_evaluation = "ENDF/BVII.1";

   allow_nat = false;
   allow_any = false;

   SetMinKinEnergy(  0*MeV );
   SetMaxKinEnergy( 20*MeV );

   lend_manager = G4LENDManager::GetInstance(); 

}
   
G4LENDCrossSection::~G4LENDCrossSection()
{

   for ( std::map< G4int , G4LENDUsedTarget* >::iterator 
         it = usedTarget_map.begin() ; it != usedTarget_map.end() ; it ++ )
   { 
      delete it->second;  
   }

}
   
void G4LENDCrossSection::BuildPhysicsTable( const G4ParticleDefinition&  )
{
   create_used_target_map();
}

void G4LENDCrossSection::DumpPhysicsTable(const G4ParticleDefinition& aP)
{

  if ( &aP != proj ) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use LEND data for particles other than neutrons!!!");  

   G4cout << G4endl;
   G4cout << "Dump Cross Sections of " << GetName() << G4endl;
   G4cout << "(Pointwise cross-section at 300 Kelvin.)" << G4endl;
   G4cout << G4endl;

   G4cout << "Target informaiton " << G4endl;

   for ( std::map< G4int , G4LENDUsedTarget* >::iterator 
         it = usedTarget_map.begin() ; it != usedTarget_map.end() ; it ++ )
   {
      G4cout 
         << "Wanted " << it->second->GetWantedEvaluation() 
         << ", Z= " << it->second->GetWantedZ() 
         << ", A= " << it->second->GetWantedA() 
         << "; Actual " << it->second->GetActualEvaluation() 
         << ", Z= " << it->second->GetActualZ() 
         << ", A= " << it->second->GetActualA() 
         << ", " << it->second->GetTarget() 
         << G4endl; 

      G4int ie = 0;

      G4GIDI_target* aTarget = it->second->GetTarget();
      G4double aT = 300;
      for ( ie = 0 ; ie < 130 ; ie++ )
      {
         G4double ke = 1.0e-5 * G4Pow::GetInstance()->powA ( 10.0 , ie/10.0 ) *eV;

         if ( ke < 20*MeV )
         {
            G4cout << "  " << GetName() << ", cross section at " << ke/eV << " [eV] = " << getLENDCrossSection ( aTarget , ke , aT )/barn << " [barn] " << G4endl;
         }
      }
      G4cout << G4endl;

   }

}


/*
//110810
//G4double G4LENDCrossSection::GetCrossSection(const G4DynamicParticle* aP , const G4Element* anElement , G4double aT)
G4double G4LENDCrossSection::GetCrossSection(const G4DynamicParticle* aP , int iZ , const G4Material* aMat)
{

//110810
   G4double aT = aMat->GetTemperature();
   G4Element* anElement = lend_manager->GetNistElementBuilder()->FindOrBuildElement( iZ );

   G4double ke = aP->GetKineticEnergy();
   G4double XS = 0.0;

   G4int numberOfIsotope = anElement->GetNumberOfIsotopes(); 

   if ( numberOfIsotope > 0 )
   {
      // User Defined Abundances   
      for ( G4int i_iso = 0 ; i_iso < numberOfIsotope ; i_iso++ )
      {

         G4int iZ = anElement->GetIsotope( i_iso )->GetZ();
         G4int iA = anElement->GetIsotope( i_iso )->GetN();
         G4double ratio = anElement->GetRelativeAbundanceVector()[i_iso];

         G4GIDI_target* aTarget = usedTarget_map.find( lend_manager->GetNucleusEncoding( iZ , iA ) )->second->GetTarget();
         XS += ratio*getLENDCrossSection ( aTarget , ke , aT );

      }
   }
   else
   {
      // Natural Abundances   
      G4NistElementBuilder* nistElementBuild = lend_manager->GetNistElementBuilder();
      G4int iZ = int ( anElement->GetZ() );
      G4int numberOfNistIso = nistElementBuild->GetNumberOfNistIsotopes( int ( anElement->GetZ() ) ); 

       G4int Nfirst = nistElementBuild->GetNistFirstIsotopeN( iZ );
      for ( G4int i = 0 ; i < numberOfNistIso ; i++ )
      {
         G4int iA = Nfirst + i;  
         G4double ratio = nistElementBuild->GetIsotopeAbundance( iZ , iA );
         if ( ratio > 0.0 )
         {
            G4GIDI_target* aTarget = usedTarget_map.find( lend_manager->GetNucleusEncoding( iZ , iA ) )->second->GetTarget();
            XS += ratio*getLENDCrossSection ( aTarget , ke , aT );
            //G4cout << ke/eV << " "  << iZ << " " << iMass << " " << aTarget << " " << getLENDCrossSection ( aTarget , ke , aT ) << G4endl;
         }
      }
   }
 
   //G4cout << "XS= " << XS << G4endl;
   return XS;
}



//110810
//G4double G4LENDCrossSection::GetIsoCrossSection(const G4DynamicParticle* dp, const G4Isotope* isotope, G4double aT )
G4double G4LENDCrossSection::GetIsoCrossSection(const G4DynamicParticle* dp, const G4Isotope* isotope, const G4Material* aMat)
{

//110810
   G4double aT = aMat->GetTemperature();

   G4double ke = dp->GetKineticEnergy();

   G4int iZ = isotope->GetZ();
   G4int iA = isotope->GetN();

   G4GIDI_target* aTarget = usedTarget_map.find( lend_manager->GetNucleusEncoding( iZ , iA ) )->second->GetTarget();

   return getLENDCrossSection ( aTarget , ke , aT );

}



//110810
//G4double G4LENDCrossSection::GetZandACrossSection(const G4DynamicParticle* dp, G4int iZ, G4int iA, G4double aT)
G4double G4LENDCrossSection::GetZandACrossSection(const G4DynamicParticle* dp, G4int iZ, G4int iA, const G4Material* aMat)
{

//110810
   G4double aT = aMat->GetTemperature();

   G4double ke = dp->GetKineticEnergy();

   G4GIDI_target* aTarget = usedTarget_map.find( lend_manager->GetNucleusEncoding( iZ , iA ) )->second->GetTarget();

   return getLENDCrossSection ( aTarget , ke , aT );

}
*/



void G4LENDCrossSection::recreate_used_target_map()
{
   for ( std::map< G4int , G4LENDUsedTarget* >::iterator 
         it = usedTarget_map.begin() ; it != usedTarget_map.end() ; it ++ )
   { 
      delete it->second;  
   }
   usedTarget_map.clear();

   create_used_target_map();
}



void G4LENDCrossSection::create_used_target_map()
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

            //G4LENDUsedTarget* aTarget = new G4LENDUsedTarget ( G4Neutron::Neutron() , default_evaluation , iZ , iA );  
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

   G4cout << "Dump UsedTarget for " << GetName() << G4endl;
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

                                                           // elow          ehigh       xs_elow      xs_ehigh      ke (ke < elow)
G4double G4LENDCrossSection::GetUltraLowEnergyExtrapolatedXS( G4double x1, G4double x2, G4double y1, G4double y2 , G4double ke )
{
   //XS propotinal to 1/v at low energy -> 1/root(E) 
   //XS = a * 1/root(E) + b  
   G4double a = ( y2 - y1 ) / ( 1/std::sqrt(x2) - 1/std::sqrt(x1) );
   G4double b = y1 - a * 1/std::sqrt(x1);
   G4double result = a * 1/std::sqrt(ke) + b;
   return result;
}
