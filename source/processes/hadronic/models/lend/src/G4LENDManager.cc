// Class Description
// Manager of LEND (Low Energy Nuclear Data) target (nucleus) 
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

#include "G4LENDManager.hh"
#include "G4HadronicException.hh"

#include "G4Neutron.hh"

G4LENDManager* G4LENDManager::endl_manager = NULL;

G4LENDManager::G4LENDManager()
{

   if(!getenv("G4LENDDATA")) 
      throw G4HadronicException(__FILE__, __LINE__, " Please setenv G4LENDDATA to point to the LEND files." );

   //G4String xmcf ( "/afs/slac.stanford.edu/package/geant4/vol49/LEND/xmcf/xmcf.map" );
   G4String xmcf = getenv("G4LENDDATA");

// for neutron

   GIDI4GEANT* axLEND = new GIDI4GEANT( 1 , xmcf );

   if ( proj_endl_map.find ( G4Neutron::Neutron() ) == proj_endl_map.end() )
   {
      proj_endl_map.insert ( std::pair < G4ParticleDefinition* , GIDI4GEANT* > ( G4Neutron::Neutron() , axLEND ) ); 
   }
      
   v_endl_target.clear();

   ionTable = new G4IonTable();
   nistElementBuilder = new G4NistElementBuilder( 0 );
}
   


G4LENDManager::~G4LENDManager()
{
   
// deleting target
   for ( std::vector < endl_target >::iterator 
         it = v_endl_target.begin() ; it != v_endl_target.end() ; it++ )
   {
        (*it).endl->freeTarget( it->target ); 
   }

// deleting endl
   for ( std::map < G4ParticleDefinition* , GIDI4GEANT* >::iterator 
         it = proj_endl_map.begin() ; it != proj_endl_map.end() ; it++ )
   {
      delete it->second;
   }

   delete ionTable;
   delete nistElementBuilder;

}



GIDI4GEANT_target* G4LENDManager::GetLENDTarget( G4ParticleDefinition* proj , G4String evaluation , G4int iZ , G4int iA , G4int iM )
{

   GIDI4GEANT_target* anLENDTarget = NULL;

   G4int iTarg = ionTable->GetNucleusEncoding( iZ , iA );
                                                     // G4double E=0.0, G4int J=0);

   for ( std::vector < endl_target >::iterator 
         it = v_endl_target.begin() ; it != v_endl_target.end() ; it++ )
   {
//    find the target
      if ( it->proj == proj && it->target_code == iTarg && it->evaluation == evaluation ) 
      {
         return it->target;
      }
   }

   if ( proj_endl_map.find ( proj ) == proj_endl_map.end() ) 
   {
      G4cout << proj->GetParticleName() << " is not supported by this LEND." << G4endl;
      return anLENDTarget; // return NULL 
   }

   GIDI4GEANT* xendl = proj_endl_map.find ( proj ) -> second; 

   if ( xendl->isThisDataAvailable( evaluation, iZ, iA , iM ) )
   {

      G4cout << evaluation << " for " << ionTable->GetIonName( iZ , iA , 0 ) << " is exist in this LEND." << G4endl;

      anLENDTarget = xendl->readTarget( evaluation , iZ , iA , iM );

      endl_target new_target; 
      new_target.endl = xendl; 
      new_target.target = anLENDTarget; 
      new_target.proj = proj;
      new_target.evaluation = evaluation;
      new_target.target_code = iTarg;
      
      v_endl_target.push_back( new_target );

//    found EXACT
      return anLENDTarget;

   }
   else 
   {
//    NO EXACT DATA (Evaluatino & Z,A,M)
                                                            // This is for ground state
      G4cout << evaluation << " for " << ionTable->GetIonName( iZ , iA , 0 ) << " is not exist in this LEND." << G4endl;
      std::vector< std::string >* available =  xendl->getNamesOfAvailableLibraries( iZ, iA , iM );
      if ( available->size() > 0 )
      {
//       EXACT Z,A,M but Evaluation is different 
         G4cout << " However you can use following evaluation(s) for the target. " << G4endl;
         std::vector< std::string >::iterator its;
         for ( its = available->begin() ; its != available->end() ; its++ ) 
            G4cout << *its << G4endl;

         G4cout << G4endl;
      }
//      
//    checking natural abundance data for Z
//
      else if ( xendl->isThisDataAvailable( evaluation, iZ, 0 , iM ) )
      {
//       EXACT natural abundance data for the evaluation 
         G4cout << " However you can use natural abundance data for the target. " << G4endl;
      }
      else
      {
         std::vector< std::string >* available_nat =  xendl->getNamesOfAvailableLibraries( iZ, 0 , iM );
//
         if ( available_nat->size() > 0 )
         {
//          EXACT natural abundance data for Z but differnet evaluation
            G4cout << " However you can use following evaluation(s) for natural abundace of the target. " << G4endl;
            std::vector< std::string >::iterator its;
            for ( its = available_nat->begin() ; its != available_nat->end() ; its++ ) 
               G4cout << *its << G4endl;
            G4cout << G4endl;
         }
      }

//    return NULL if exact data is not available               
      return anLENDTarget; // return NULL   
   }
   
   return anLENDTarget; 
}



std::vector< G4String > G4LENDManager::IsLENDTargetAvailable ( G4ParticleDefinition* proj , G4int iZ , G4int iA , G4int iM )
{

   std::vector< G4String > answer;
   if ( proj_endl_map.find ( proj ) == proj_endl_map.end() ) 
   {
      G4cout << proj->GetParticleName() << " is not supported by this LEND." << G4endl;
      return answer; // return NULL 
   }

   GIDI4GEANT* xendl = proj_endl_map.find ( proj ) -> second; 
   std::vector< std::string >* available =  xendl->getNamesOfAvailableLibraries( iZ, iA , iM );

   if ( available->size() > 0 )
   {
      std::vector< std::string >::iterator its;
      for ( its = available->begin() ; its != available->end() ; its++ ) 
         answer.push_back ( *its );
   }

   return answer;
}
