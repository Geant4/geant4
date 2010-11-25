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

using namespace std;
#include "G4LENDCrossSection.hh"

#include "G4ElementTable.hh"
#include "G4HadronicException.hh"

G4bool G4LENDCrossSection::IsApplicable(const G4DynamicParticle*aP, const G4Element*)
{

   G4bool result = true;
   G4double eKin = aP->GetKineticEnergy();
   if( eKin > 20*MeV || aP->GetDefinition() != proj ) result = false;
   return result;

}

G4LENDCrossSection::G4LENDCrossSection()
{

   default_evaluation = "endl99";
   allow_nat = false;
   allow_any = false;

   endl_manager = G4LENDManager::GetInstance(); 

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
   ;
}

void G4LENDCrossSection::DumpPhysicsTable(const G4ParticleDefinition& aP)
{

  if ( &aP != proj ) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use LEND data for particles other than neutrons!!!");  

   G4cout << G4endl;
   G4cout << "Dump Cross Sections of " << name << G4endl;
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

      GIDI4GEANT_target* aTarget = it->second->GetTarget();
      G4double aT = 300;
      for ( ie = 0 ; ie < 130 ; ie++ )
      {
         G4double ke = 1.0e-5 * std::pow ( 10.0 , ie/10.0 ) *eV;

         if ( ke < 20*MeV )
         {
            G4cout << "  "<< name << ", cross section at " << ke/eV << " [eV] = " << getLENDCrossSection ( aTarget , ke , aT )/barn << " [barn] " << G4endl;
         }
      }
      G4cout << G4endl;

   }

}



G4double G4LENDCrossSection::GetCrossSection(const G4DynamicParticle* aP , const G4Element* anElement , G4double aT)
{

// G4cout << "G4LENDCrossSection::GetCrossSection(const G4DynamicParticle* aP, const G4Element*anE, G4double aT) " << G4endl;

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

         GIDI4GEANT_target* aTarget = usedTarget_map.find( endl_manager->GetNucleusEncoding( iZ , iA ) )->second->GetTarget();
         XS += ratio*getLENDCrossSection ( aTarget , ke , aT );

      }
   }
   else
   {
      // Natural Abundances   
      G4NistElementBuilder* nistElementBuild = endl_manager->GetNistElementBuilder();
      G4int iZ = int ( anElement->GetZ() );
      G4int numberOfNistIso = nistElementBuild->GetNumberOfNistIsotopes( int ( anElement->GetZ() ) ); 

      for ( G4int i = 0 ; i < numberOfNistIso ; i++ )
      {

         if ( nistElementBuild->GetIsotopeAbundance( iZ , nistElementBuild->GetNistFirstIsotopeN( iZ ) + i ) > 0 )
         {
            G4int iMass = nistElementBuild->GetNistFirstIsotopeN( iZ ) + i;  

            G4double ratio = nistElementBuild->GetIsotopeAbundance ( iZ , iMass );
            GIDI4GEANT_target* aTarget = usedTarget_map.find( endl_manager->GetNucleusEncoding( iZ , iMass ) )->second->GetTarget();
            XS += ratio*getLENDCrossSection ( aTarget , ke , aT );

         }

      }
   }
 
   return XS;
}



G4double G4LENDCrossSection::GetIsoCrossSection(const G4DynamicParticle* dp, const G4Isotope* isotope, G4double aT )
{

   G4double ke = dp->GetKineticEnergy();

   G4int iZ = isotope->GetZ();
   G4int iA = isotope->GetN();

   GIDI4GEANT_target* aTarget = usedTarget_map.find( endl_manager->GetNucleusEncoding( iZ , iA ) )->second->GetTarget();

   return getLENDCrossSection ( aTarget , ke , aT );

}



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

            //G4LENDUsedTarget* aTarget = new G4LENDUsedTarget ( G4Neutron::Neutron() , default_evaluation , iZ , iA );  
            G4LENDUsedTarget* aTarget = new G4LENDUsedTarget ( proj , default_evaluation , iZ , iA );  
            if ( allow_nat == true ) aTarget->AllowNat();
            if ( allow_any == true ) aTarget->AllowAny();
            usedTarget_map.insert( std::pair< G4int , G4LENDUsedTarget* > ( endl_manager->GetNucleusEncoding( iZ , iA ) , aTarget ) );
         }
      }
      else
      {
      // Natural Abundances   
         G4NistElementBuilder* nistElementBuild = endl_manager->GetNistElementBuilder();
         G4int iZ = int ( anElement->GetZ() );
         //G4cout << nistElementBuild->GetNumberOfNistIsotopes( int ( anElement->GetZ() ) ) << G4endl;
         G4int numberOfNistIso = nistElementBuild->GetNumberOfNistIsotopes( int ( anElement->GetZ() ) ); 

         for ( G4int i = 0 ; i < numberOfNistIso ; i++ )
         {
            //G4cout << nistElementBuild->GetIsotopeAbundance( iZ , nistElementBuild->GetNistFirstIsotopeN( iZ ) + i ) << G4endl;
            if ( nistElementBuild->GetIsotopeAbundance( iZ , nistElementBuild->GetNistFirstIsotopeN( iZ ) + i ) > 0 )
            {
               G4int iMass = nistElementBuild->GetNistFirstIsotopeN( iZ ) + i;  
               //G4cout << iZ << " " << nistElementBuild->GetNistFirstIsotopeN( iZ ) + i << " " << nistElementBuild->GetIsotopeAbundance ( iZ , iMass ) << G4endl;  

               G4LENDUsedTarget* aTarget = new G4LENDUsedTarget ( proj , default_evaluation , iZ , iMass );  
               if ( allow_nat == true ) aTarget->AllowNat();
               if ( allow_any == true ) aTarget->AllowAny();
               usedTarget_map.insert( std::pair< G4int , G4LENDUsedTarget* > ( endl_manager->GetNucleusEncoding( iZ , iMass ) , aTarget ) );

            }

         }
      }
   }

   G4cout << "Dump UsedTarget for " << name << G4endl;
   for ( std::map< G4int , G4LENDUsedTarget* >::iterator 
         it = usedTarget_map.begin() ; it != usedTarget_map.end() ; it ++ )
   {
      G4cout 
         << " " << it->second->GetWantedEvaluation() 
         << " " << it->second->GetWantedZ() 
         << " " << it->second->GetWantedA() 
         << " " << it->second->GetActualEvaluation() 
         << " " << it->second->GetActualZ() 
         << " " << it->second->GetActualA() 
         << " " << it->second->GetTarget() 
         << G4endl; 
   } 

}
