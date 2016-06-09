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

#include "G4NeutronHPThermalScatteringData.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
//#include "G4NeutronHPData.hh"



G4NeutronHPThermalScatteringData::G4NeutronHPThermalScatteringData()
{
// Upper limit of neutron energy 
   emax = 4*eV;

   indexOfThermalElement.clear(); 

   names = new G4NeutronHPThermalScatteringNames();

   BuildPhysicsTable( *G4Neutron::Neutron() );
}



G4NeutronHPThermalScatteringData::~G4NeutronHPThermalScatteringData()
{

   clearCurrentXSData();

   delete names;
}


void G4NeutronHPThermalScatteringData::clearCurrentXSData()
{
   std::map< G4int , std::map< G4double , G4NeutronHPVector* >* >::iterator it;
   std::map< G4double , G4NeutronHPVector* >::iterator itt;

   for ( it = coherent.begin() ; it != coherent.end() ; it++ )
   {
      if ( it->second != NULL )
      {
         for ( itt = it->second->begin() ; itt != it->second->end() ; itt++ )
         {
            delete itt->second;
         }
      }
      delete it->second;
   }

   for ( it = incoherent.begin() ; it != incoherent.end() ; it++ )
   {
      if ( it->second != NULL )
      { 
         for ( itt = it->second->begin() ; itt != it->second->end() ; itt++ )
         {
            delete itt->second;
         }
      }
      delete it->second;
   }

   for ( it = inelastic.begin() ; it != inelastic.end() ; it++ )
   {
      if ( it->second != NULL )
      {
         for ( itt = it->second->begin() ; itt != it->second->end() ; itt++ )
         {
            delete itt->second;
         }
      }
      delete it->second; 
   }

   coherent.clear();
   incoherent.clear();
   inelastic.clear();

}



G4bool G4NeutronHPThermalScatteringData::IsApplicable( const G4DynamicParticle* aP , const G4Element* anEle )
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
         std::vector < G4int >::iterator it; 
         for ( it = indexOfThermalElement.begin() ; it != indexOfThermalElement.end() ; it++ )
         {
             if ( ie == *it ) return true;
         }
      }
   }

/*
   if ( names->IsThisThermalElement ( anEle->GetName() ) )
   {
      // Check energy and projectile species 
      G4double eKin = aP->GetKineticEnergy();
      if ( eKin < emax && aP->GetDefinition() == G4Neutron::Neutron() ) result = true; 
   }
*/
   return result;
}


void G4NeutronHPThermalScatteringData::BuildPhysicsTable(const G4ParticleDefinition& aP)
{

   if ( &aP != G4Neutron::Neutron() ) 
      throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  

   indexOfThermalElement.clear(); 

   clearCurrentXSData();

   static const G4ElementTable* theElementTable = G4Element::GetElementTable();
   size_t numberOfElements = G4Element::GetNumberOfElements();
   size_t numberOfThermalElements = 0; 
   for ( size_t i = 0 ; i < numberOfElements ; i++ )
   {
      if ( names->IsThisThermalElement ( (*theElementTable)[i]->GetName() ) )
      {
         indexOfThermalElement.push_back( i ); 
         numberOfThermalElements++;
      }
   }

   // Read Cross Section Data files

   G4String dirName;
   if ( !getenv( "G4NEUTRONHPDATA" ) ) 
      throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
   G4String baseName = getenv( "G4NEUTRONHPDATA" );

   dirName = baseName + "/ThermalScattering";

   G4String ndl_filename;
   G4String name;

   for ( size_t i = 0 ; i < numberOfThermalElements ; i++ )
   {
      ndl_filename = names->GetTS_NDL_Name( (*theElementTable)[ indexOfThermalElement[ i ] ]->GetName() ); 

      // Coherent
      name = dirName + "/Coherent/CrossSection/" + ndl_filename; 
      std::map< G4double , G4NeutronHPVector* >*  coh_amapTemp_EnergyCross = readData( name );
      coherent.insert ( std::pair < G4int , std::map< G4double , G4NeutronHPVector* >* > ( indexOfThermalElement[ i ] , coh_amapTemp_EnergyCross ) );

      // Incoherent
      name = dirName + "/Incoherent/CrossSection/" + ndl_filename; 
      std::map< G4double , G4NeutronHPVector* >*  incoh_amapTemp_EnergyCross = readData( name );
      incoherent.insert ( std::pair < G4int , std::map< G4double , G4NeutronHPVector* >* > ( indexOfThermalElement[ i ] , incoh_amapTemp_EnergyCross ) );

      // Inelastic
      name = dirName + "/Inelastic/CrossSection/" + ndl_filename; 
      std::map< G4double , G4NeutronHPVector* >*  inela_amapTemp_EnergyCross = readData( name );
      inelastic.insert ( std::pair < G4int , std::map< G4double , G4NeutronHPVector* >* > ( indexOfThermalElement[ i ] , inela_amapTemp_EnergyCross ) );
   }

}



std::map< G4double , G4NeutronHPVector* >* G4NeutronHPThermalScatteringData::readData ( G4String name ) 
{

   std::map< G4double , G4NeutronHPVector* >*  aData = new std::map< G4double , G4NeutronHPVector* >; 
   
   std::ifstream theChannel( name.c_str() );

   //G4cout << "G4NeutronHPThermalScatteringData " << name << G4endl;

   G4int dummy; 
   while ( theChannel >> dummy )   // MF
   {
      theChannel >> dummy;   // MT
      G4double temp; 
      theChannel >> temp;   
      G4NeutronHPVector* anEnergyCross = new G4NeutronHPVector;
      G4int nData;
      theChannel >> nData;
      anEnergyCross->Init ( theChannel , nData , eV , barn );
      aData->insert ( std::pair < G4double , G4NeutronHPVector* > ( temp , anEnergyCross ) );
   }
   theChannel.close();

   return aData;

} 



void G4NeutronHPThermalScatteringData::DumpPhysicsTable( const G4ParticleDefinition& aP )
{
   if( &aP != G4Neutron::Neutron() ) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  
//  G4cout << "G4NeutronHPThermalScatteringData::DumpPhysicsTable still to be implemented"<<G4endl;
}

//#include "G4Nucleus.hh"
//#include "G4NucleiPropertiesTable.hh"
//#include "G4Neutron.hh"
//#include "G4Electron.hh"



G4double G4NeutronHPThermalScatteringData::GetCrossSection( const G4DynamicParticle* aP , const G4Element*anE , G4double aT )
{
   G4double result = 0;
   
   G4int iele = anE->GetIndex();

   G4double Xcoh = GetX ( aP , aT , coherent.find(iele)->second );
   G4double Xincoh = GetX ( aP , aT , incoherent.find(iele)->second );
   G4double Xinela = GetX ( aP , aT , inelastic.find(iele)->second );

   result = Xcoh + Xincoh + Xinela;

   //G4cout << "G4NeutronHPThermalScatteringData::GetCrossSection  Tot= " << result/barn << " Coherent= " << Xcoh/barn << " Incoherent= " << Xincoh/barn << " Inelastic= " << Xinela/barn << G4endl;

   return result;
}



G4double G4NeutronHPThermalScatteringData::GetInelasticCrossSection( const G4DynamicParticle* aP , const G4Element*anE , G4double aT )
{
   G4double result = 0;
   G4int iele = anE->GetIndex();
   result = GetX ( aP , aT , inelastic.find(iele)->second );
   return result;
}



G4double G4NeutronHPThermalScatteringData::GetCoherentCrossSection( const G4DynamicParticle* aP , const G4Element*anE , G4double aT )
{
   G4double result = 0;
   G4int iele = anE->GetIndex();
   result = GetX ( aP , aT , coherent.find(iele)->second );
   return result;
}



G4double G4NeutronHPThermalScatteringData::GetIncoherentCrossSection( const G4DynamicParticle* aP , const G4Element*anE , G4double aT )
{
   G4double result = 0;
   G4int iele = anE->GetIndex();
   result = GetX ( aP , aT , incoherent.find(iele)->second );
   return result;
}




G4double G4NeutronHPThermalScatteringData::GetX ( const G4DynamicParticle* aP, G4double aT , std::map < G4double , G4NeutronHPVector* >* amapTemp_EnergyCross )
{
   G4double result = 0;
   if ( amapTemp_EnergyCross->size() == 0 ) return result;

   std::map< G4double , G4NeutronHPVector* >::iterator it; 
   for ( it = amapTemp_EnergyCross->begin() ; it != amapTemp_EnergyCross->end() ; it++ )
   {
       if ( aT < it->first ) break;
   } 
   if ( it == amapTemp_EnergyCross->begin() ) it++;  // lower than first
      else if ( it == amapTemp_EnergyCross->end() ) it--;  // upper than last

   G4double eKinetic = aP->GetKineticEnergy();

   G4double TH = it->first;
   G4double XH = it->second->GetXsec ( eKinetic ); 

   //G4cout << "G4NeutronHPThermalScatteringData::GetX TH " << TH << " E " << eKinetic <<  " XH " << XH << G4endl;

   it--;
   G4double TL = it->first;
   G4double XL = it->second->GetXsec ( eKinetic ); 

   //G4cout << "G4NeutronHPThermalScatteringData::GetX TL " << TL << " E " << eKinetic <<  " XL " << XL << G4endl;

   if ( TH == TL )  
      throw G4HadronicException(__FILE__, __LINE__, "Thermal Scattering Data Error!");  

   G4double T = aT;
   G4double X = ( XH - XL ) / ( TH - TL ) * ( T - TL ) + XL;
   result = X;
  
   return result;
}
