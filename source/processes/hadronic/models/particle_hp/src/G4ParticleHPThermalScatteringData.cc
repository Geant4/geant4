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
// Thermal Neutron Scattering
// Koi, Tatsumi (SCCS/SLAC)
//
// Class Description
// Cross Sections for a high precision (based on evaluated data
// libraries) description of themal neutron scattering below 4 eV;
// Based on Thermal neutron scattering files
// from the evaluated nuclear data files ENDF/B-VI, Release2
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with
// the corresponding process.
// Class Description - End

// 15-Nov-06 First implementation is done by T. Koi (SLAC/SCCS)
// 070625 implement clearCurrentXSData to fix memory leaking by T. Koi
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//

#include <list>
#include <algorithm>

#include "G4ParticleHPThermalScatteringData.hh"
#include "G4ParticleHPManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"

#include "G4Threading.hh"

G4ParticleHPThermalScatteringData::G4ParticleHPThermalScatteringData()
:G4VCrossSectionDataSet("NeutronHPThermalScatteringData")
,coherent(nullptr)
,incoherent(nullptr)
,inelastic(nullptr)
{
   // Upper limit of neutron energy 
   emax = 4*eV;
   SetMinKinEnergy( 0*MeV );                                   
   SetMaxKinEnergy( emax );                                   

   ke_cache = 0.0;
   xs_cache = 0.0;
   element_cache = nullptr;
   material_cache = nullptr;

   indexOfThermalElement.clear(); 

   names = new G4ParticleHPThermalScatteringNames();
}

G4ParticleHPThermalScatteringData::~G4ParticleHPThermalScatteringData()
{
   clearCurrentXSData();

   delete names;
}

G4bool G4ParticleHPThermalScatteringData::IsIsoApplicable( const G4DynamicParticle* dp , 
                                                G4int /*Z*/ , G4int /*A*/ ,
                                                const G4Element* element ,
                                                const G4Material* material )
{
   G4double eKin = dp->GetKineticEnergy();
   if ( eKin > 4.0*eV //GetMaxKinEnergy() 
     || eKin < 0 //GetMinKinEnergy() 
     || dp->GetDefinition() != G4Neutron::Neutron() ) return false;                                   

   if ( dic.find( std::pair < const G4Material* , const G4Element* > ( (G4Material*)NULL , element ) ) != dic.end() 
     || dic.find( std::pair < const G4Material* , const G4Element* > ( material , element ) ) != dic.end() ) return true;

   return false;
}

G4double G4ParticleHPThermalScatteringData::GetIsoCrossSection( const G4DynamicParticle* dp ,
                                   G4int /*Z*/ , G4int /*A*/ ,
                                   const G4Isotope* /*iso*/  ,
                                   const G4Element* element ,
                                   const G4Material* material )
{
   ke_cache = dp->GetKineticEnergy();
   element_cache = element;
   material_cache = material;
   G4double xs = GetCrossSection( dp , element , material );
   xs_cache = xs;
   return xs;
}

void G4ParticleHPThermalScatteringData::clearCurrentXSData()
{
   if ( coherent != nullptr )
   {
     for (auto it=coherent->cbegin() ; it!=coherent->cend(); ++it)
     {
       if ( it->second != nullptr )
       {
         for (auto itt=it->second->cbegin(); itt!=it->second->cend(); ++itt)
         {
           delete itt->second;
         }
       }
       delete it->second;
     }
     coherent->clear();
   }

   if ( incoherent != nullptr )
   {
     for (auto it=incoherent->cbegin(); it!=incoherent->cend() ; ++it)
     {
       if ( it->second != nullptr )
       {
         for (auto itt=it->second->cbegin(); itt!=it->second->cend(); ++itt)
         {
           delete itt->second;
         }
       }
       delete it->second;
     }
     incoherent->clear();
   }

   if ( inelastic != nullptr )
   {
     for (auto it=inelastic->cbegin(); it!=inelastic->cend(); ++it)
     {
       if ( it->second != nullptr )
       {
         for (auto itt=it->second->cbegin(); itt!=it->second->cend(); ++itt)
         {
           delete itt->second;
         }
       }
       delete it->second;
     }
     inelastic->clear();
   }

}


G4bool G4ParticleHPThermalScatteringData::IsApplicable( const G4DynamicParticle* aP , const G4Element* anEle )
{
   G4bool result = false;

   G4double eKin = aP->GetKineticEnergy();
   // Check energy 
   if ( eKin < emax )
   {
      // Check Particle Species
      if ( aP->GetDefinition() == G4Neutron::Neutron() ) 
      {
        // anEle is one of Thermal elements 
         G4int ie = (G4int) anEle->GetIndex();
         for (auto it = indexOfThermalElement.cbegin();
                   it != indexOfThermalElement.cend() ; ++it)
         {
             if ( ie == *it ) return true;
         }
      }
   }

   return result;
}


void G4ParticleHPThermalScatteringData::BuildPhysicsTable(const G4ParticleDefinition& aP)
{

   if ( &aP != G4Neutron::Neutron() ) 
      throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  

   //std::map < std::pair < G4Material* , const G4Element* > , G4int > dic;   
   //
   dic.clear();   
   if ( G4Threading::IsMasterThread() ) clearCurrentXSData();

   std::map < G4String , G4int > co_dic;   

   //Searching Nist Materials
   static G4ThreadLocal G4MaterialTable* theMaterialTable  = nullptr;
   if (!theMaterialTable) theMaterialTable= G4Material::GetMaterialTable();
   std::size_t numberOfMaterials = G4Material::GetNumberOfMaterials();
   for ( std::size_t i = 0 ; i < numberOfMaterials ; ++i )
   {
      G4Material* material = (*theMaterialTable)[i];
      G4int numberOfElements = (G4int)material->GetNumberOfElements();
      for ( G4int j = 0 ; j < numberOfElements ; ++j )
      {
         const G4Element* element = material->GetElement(j);
         if ( names->IsThisThermalElement ( material->GetName() , element->GetName() ) )
         {                                    
            G4int ts_ID_of_this_geometry; 
            G4String ts_ndl_name = names->GetTS_NDL_Name( material->GetName() , element->GetName() ); 
            if ( co_dic.find ( ts_ndl_name ) != co_dic.cend() )
            {
               ts_ID_of_this_geometry = co_dic.find ( ts_ndl_name ) -> second;
            }
            else
            {
               ts_ID_of_this_geometry = (G4int)co_dic.size();
               co_dic.insert ( std::pair< G4String , G4int >( ts_ndl_name , ts_ID_of_this_geometry ) );
            }

            dic.insert( std::pair < std::pair < G4Material* , const G4Element* > , G4int > ( std::pair < G4Material* , const G4Element* > ( material , element ) , ts_ID_of_this_geometry ) );
         }
      }
   }

   //Searching TS Elements 
   static G4ThreadLocal G4ElementTable* theElementTable  = nullptr;
   if (!theElementTable) theElementTable= G4Element::GetElementTable();
   std::size_t numberOfElements = G4Element::GetNumberOfElements();

   for ( std::size_t i = 0 ; i < numberOfElements ; ++i )
   {
      const G4Element* element = (*theElementTable)[i];
      if ( names->IsThisThermalElement ( element->GetName() ) )
      {
         if ( names->IsThisThermalElement ( element->GetName() ) )
         {                                    
            G4int ts_ID_of_this_geometry; 
            G4String ts_ndl_name = names->GetTS_NDL_Name( element->GetName() ); 
            if ( co_dic.find ( ts_ndl_name ) != co_dic.cend() )
            {
               ts_ID_of_this_geometry = co_dic.find ( ts_ndl_name ) -> second;
            }
            else
            {
               ts_ID_of_this_geometry = (G4int)co_dic.size();
               co_dic.insert ( std::pair< G4String , G4int >( ts_ndl_name , ts_ID_of_this_geometry ) );
            }

            dic.insert( std::pair < std::pair < const G4Material* , const G4Element* > , G4int > ( std::pair < const G4Material* , const G4Element* > ( (G4Material*)NULL , element ) ,  ts_ID_of_this_geometry ) );
         }
      }
   }

   G4cout << G4endl;
   G4cout << "Neutron HP Thermal Scattering Data: Following material-element pairs and/or elements are registered." << G4endl;
   for ( auto it = dic.cbegin() ; it != dic.cend() ; ++it )   
   {
      if ( it->first.first != nullptr ) 
      {
         G4cout << "Material " << it->first.first->GetName() << " - Element "
                << it->first.second->GetName()
                << ",  internal thermal scattering id " << it->second << G4endl;
      }
      else
      {
         G4cout << "Element " << it->first.second->GetName()
                << ",  internal thermal scattering id " << it->second << G4endl;
      }
   }
   G4cout << G4endl;

   G4ParticleHPManager* hpmanager = G4ParticleHPManager::GetInstance();

   coherent = hpmanager->GetThermalScatteringCoherentCrossSections();
   incoherent = hpmanager->GetThermalScatteringIncoherentCrossSections();
   inelastic = hpmanager->GetThermalScatteringInelasticCrossSections();

   if ( G4Threading::IsMasterThread() )
   {
      if ( coherent == nullptr )
        coherent = new std::map< G4int , std::map< G4double , G4ParticleHPVector* >* >;
      if ( incoherent == nullptr )
        incoherent = new std::map< G4int , std::map< G4double , G4ParticleHPVector* >* >;
      if ( inelastic == nullptr )
        inelastic = new std::map< G4int , std::map< G4double , G4ParticleHPVector* >* >;

      // Read Cross Section Data files

      G4String dirName;
      if ( !G4FindDataDir( "G4NEUTRONHPDATA" ) )
         throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
      G4String baseName = G4FindDataDir( "G4NEUTRONHPDATA" );

      dirName = baseName + "/ThermalScattering";

      G4String ndl_filename;
      G4String full_name;

      for ( auto it = co_dic.cbegin() ; it != co_dic.cend() ; ++it )  
      {
         ndl_filename = it->first;
         G4int ts_ID = it->second;

         // Coherent
         full_name = dirName + "/Coherent/CrossSection/" + ndl_filename; 
         auto  coh_amapTemp_EnergyCross = readData( full_name );
         coherent->insert ( std::pair < G4int , std::map< G4double , G4ParticleHPVector* >* > ( ts_ID , coh_amapTemp_EnergyCross ) );

         // Incoherent
         full_name = dirName + "/Incoherent/CrossSection/" + ndl_filename; 
         auto  incoh_amapTemp_EnergyCross = readData( full_name );
         incoherent->insert ( std::pair < G4int , std::map< G4double , G4ParticleHPVector* >* > ( ts_ID , incoh_amapTemp_EnergyCross ) );

         // Inelastic
         full_name = dirName + "/Inelastic/CrossSection/" + ndl_filename; 
         auto  inela_amapTemp_EnergyCross = readData( full_name );
         inelastic->insert ( std::pair < G4int , std::map< G4double , G4ParticleHPVector* >* > ( ts_ID , inela_amapTemp_EnergyCross ) );
      }
      hpmanager->RegisterThermalScatteringCoherentCrossSections( coherent );
      hpmanager->RegisterThermalScatteringIncoherentCrossSections( incoherent );
      hpmanager->RegisterThermalScatteringInelasticCrossSections( inelastic );
   } 
}


std::map< G4double , G4ParticleHPVector* >*
G4ParticleHPThermalScatteringData::readData ( G4String full_name ) 
{
   auto  aData = new std::map< G4double , G4ParticleHPVector* >; 
   
   std::istringstream theChannel;
   G4ParticleHPManager::GetInstance()->GetDataStream(full_name,theChannel);

   G4int dummy; 
   while ( theChannel >> dummy )   // MF // Loop checking, 11.05.2015, T. Koi
   {
      theChannel >> dummy;   // MT
      G4double temp; 
      theChannel >> temp;   
      G4ParticleHPVector* anEnergyCross = new G4ParticleHPVector;
      G4int nData;
      theChannel >> nData;
      anEnergyCross->Init ( theChannel , nData , eV , barn );
      aData->insert ( std::pair < G4double , G4ParticleHPVector* > ( temp , anEnergyCross ) );
   }

   return aData;
} 


void G4ParticleHPThermalScatteringData::DumpPhysicsTable( const G4ParticleDefinition& aP )
{
   if( &aP != G4Neutron::Neutron() ) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  
}


G4double G4ParticleHPThermalScatteringData::GetCrossSection( const G4DynamicParticle* aP , const G4Element*anE , const G4Material* aM )
{
   G4double result = 0;
   
   G4int ts_id =getTS_ID( aM , anE );

   if ( ts_id == -1 ) return result;

   G4double aT = aM->GetTemperature();

   G4double Xcoh = GetX ( aP , aT , coherent->find(ts_id)->second );
   G4double Xincoh = GetX ( aP , aT , incoherent->find(ts_id)->second );
   G4double Xinela = GetX ( aP , aT , inelastic->find(ts_id)->second );

   result = Xcoh + Xincoh + Xinela;

   return result;
}


G4double G4ParticleHPThermalScatteringData::GetInelasticCrossSection( const G4DynamicParticle* aP , const G4Element*anE , const G4Material* aM )
{
   G4double result = 0;
   G4int ts_id = getTS_ID( aM , anE );
   G4double aT = aM->GetTemperature();
   result = GetX ( aP , aT , inelastic->find( ts_id )->second );
   return result;
}

G4double G4ParticleHPThermalScatteringData::GetCoherentCrossSection( const G4DynamicParticle* aP , const G4Element*anE , const G4Material* aM )
{
   G4double result = 0;
   G4int ts_id = getTS_ID( aM , anE );
   G4double aT = aM->GetTemperature();
   result = GetX ( aP , aT , coherent->find( ts_id )->second );
   return result;
}

G4double G4ParticleHPThermalScatteringData::GetIncoherentCrossSection( const G4DynamicParticle* aP , const G4Element*anE , const G4Material* aM )
{
   G4double result = 0;
   G4int ts_id = getTS_ID( aM , anE );
   G4double aT = aM->GetTemperature();
   result = GetX ( aP , aT , incoherent->find( ts_id )->second );
   return result;
}

G4int G4ParticleHPThermalScatteringData::getTS_ID ( const G4Material* material , const G4Element* element )
{
   G4int result = -1;
   if ( dic.find( std::pair < const G4Material* , const G4Element* > ( (G4Material*)NULL , element ) ) != dic.end() ) 
      return dic.find( std::pair < const G4Material* , const G4Element* > ( (G4Material*)NULL , element ) )->second; 
   if ( dic.find( std::pair < const G4Material* , const G4Element* > ( material , element ) ) != dic.end() ) 
      return dic.find( std::pair < const G4Material* , const G4Element* > ( material , element ) )->second; 
   return result; 
}

G4double G4ParticleHPThermalScatteringData::GetX ( const G4DynamicParticle* aP, G4double aT , std::map < G4double , G4ParticleHPVector* >* amapTemp_EnergyCross )
{
   G4double result = 0;
   if ( amapTemp_EnergyCross->size() == 0 ) return result;

   G4double eKinetic = aP->GetKineticEnergy();

   if ( amapTemp_EnergyCross->size() == 1 ) { 
      if ( std::fabs ( aT - amapTemp_EnergyCross->cbegin()->first ) / amapTemp_EnergyCross->begin()->first > 0.1 ) {
         G4cout << "G4ParticleHPThermalScatteringData:: The temperature of material (" 
                << aT/kelvin << "K) is different more than 10% from temperature of thermal scattering file expected (" 
                << amapTemp_EnergyCross->begin()->first << "K). Result may not be reliable."
         << G4endl;
      }
      result = amapTemp_EnergyCross->begin()->second->GetXsec ( eKinetic ); 
      return result;
   }

   auto it = amapTemp_EnergyCross->cbegin();
   for (it=amapTemp_EnergyCross->cbegin(); it!=amapTemp_EnergyCross->cend(); ++it)
   {
      if ( aT < it->first ) break;
   } 
   if ( it == amapTemp_EnergyCross->cbegin() ) {
      ++it;  // lower than the first
   } else if ( it == amapTemp_EnergyCross->cend() ) {
      --it;  // upper than the last
   }

   G4double TH = it->first;
   G4double XH = it->second->GetXsec ( eKinetic ); 

   if ( it != amapTemp_EnergyCross->cbegin() ) --it;
   G4double TL = it->first;
   G4double XL = it->second->GetXsec ( eKinetic ); 

   if ( TH == TL )  
      throw G4HadronicException(__FILE__, __LINE__, "Thermal Scattering Data Error!");  

   G4double T = aT;
   G4double X = ( XH - XL ) / ( TH - TL ) * ( T - TL ) + XL;
   result = X;
  
   return result;
}


void G4ParticleHPThermalScatteringData::AddUserThermalScatteringFile( G4String nameG4Element , G4String filename )
{
   names->AddThermalElement( nameG4Element , filename );
}

void G4ParticleHPThermalScatteringData::CrossSectionDescription(std::ostream& outFile) const
{
    outFile << "High Precision cross data based on thermal scattering data in evaluated nuclear data libraries for neutrons below 5eV on specific materials\n" ;
}
