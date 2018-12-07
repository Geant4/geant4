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
// Manager of NetronHP
// 
// 121031 First implementation done by T. Koi (SLAC/PPA)
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include <zlib.h>
#include <fstream>

#include "G4ParticleHPManager.hh"
#include "G4ParticleHPThreadLocalManager.hh"
#include "G4ParticleHPMessenger.hh"
#include "G4HadronicException.hh"

G4ParticleHPManager* G4ParticleHPManager::instance = 0;

G4ParticleHPManager::G4ParticleHPManager()
: verboseLevel(1)
,USE_ONLY_PHOTONEVAPORATION(false)
,SKIP_MISSING_ISOTOPES(false)
,NEGLECT_DOPPLER(false)
,DO_NOT_ADJUST_FINAL_STATE(false)
,PRODUCE_FISSION_FRAGMENTS(false)
,USE_NRESP71_MODEL(false)
,theElasticCrossSections(0)
,theCaptureCrossSections(0)
,theFissionCrossSections(0)
,theElasticFSs(0)
,theCaptureFSs(0)
,theFissionFSs(0)
,theTSCoherentCrossSections(0)
,theTSIncoherentCrossSections(0)
,theTSInelasticCrossSections(0)
,theTSCoherentFinalStates(0)
,theTSIncoherentFinalStates(0)
,theTSInelasticFinalStates(0)
{
   messenger = new G4ParticleHPMessenger( this );
   if ( getenv( "G4NEUTRONHP_DO_NOT_ADJUST_FINAL_STATE" ) || getenv("G4PHP_DO_NOT_ADJUST_FINAL_STATE") ) DO_NOT_ADJUST_FINAL_STATE = true;
   if ( getenv( "G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION" ) ) USE_ONLY_PHOTONEVAPORATION = true;
   if ( getenv( "G4NEUTRONHP_NEGLECT_DOPPLER" ) || getenv("G4PHP_NEGLECT_DOPPLER") ) NEGLECT_DOPPLER = true;
   if ( getenv( "G4NEUTRONHP_SKIP_MISSING_ISOTOPES" ) ) SKIP_MISSING_ISOTOPES = true;
   if ( getenv( "G4NEUTRONHP_PRODUCE_FISSION_FRAGMENTS" ) ) PRODUCE_FISSION_FRAGMENTS = true;
   if ( getenv( "G4PHP_USE_NRESP71_MODEL" ) ) USE_NRESP71_MODEL = true;
}

G4ParticleHPManager::~G4ParticleHPManager()
{
   delete messenger;
}

G4ParticleHPManager* G4ParticleHPManager::GetInstance()
{
  static G4ParticleHPManager manager;
  if (!instance)
  {
    instance = &manager;
  }
  return instance;
}

void G4ParticleHPManager::OpenReactionWhiteBoard()
{
   G4ParticleHPThreadLocalManager::GetInstance()->OpenReactionWhiteBoard();
}

G4ParticleHPReactionWhiteBoard* G4ParticleHPManager::GetReactionWhiteBoard()
{
   return G4ParticleHPThreadLocalManager::GetInstance()->GetReactionWhiteBoard();
}

void G4ParticleHPManager::CloseReactionWhiteBoard()
{
   G4ParticleHPThreadLocalManager::GetInstance()->CloseReactionWhiteBoard();
}

void G4ParticleHPManager::GetDataStream( G4String filename , std::istringstream& iss ) 
{
   G4String* data=0;
   G4String compfilename(filename);
   compfilename += ".z";
   std::ifstream* in = new std::ifstream ( compfilename , std::ios::binary | std::ios::ate );
   if ( in->good() )
   {
      // Use the compressed file 
      G4int file_size = in->tellg();
      in->seekg( 0 , std::ios::beg );
      Bytef* compdata = new Bytef[ file_size ];

      while ( *in )
      { // Loop checking, 11.05.2015, T. Koi
         in->read( (char*)compdata , file_size );
      }

      uLongf complen = (uLongf) ( file_size*4 );
      Bytef* uncompdata = new Bytef[complen];

      while ( Z_OK != uncompress ( uncompdata , &complen , compdata , file_size ) )
      { // Loop checking, 11.05.2015, T. Koi
         delete[] uncompdata;
         complen *= 2;
         uncompdata = new Bytef[complen];
      }
      delete [] compdata;
      //                                 Now "complen" has uncomplessed size
      data = new G4String ( (char*)uncompdata , (G4long)complen );
      delete [] uncompdata;
   }
   else
   {
      // Use regular text file 
      std::ifstream thefData( filename , std::ios::in | std::ios::ate );
      if ( thefData.good() )
      {
         G4int file_size = thefData.tellg();
         thefData.seekg( 0 , std::ios::beg );
         char* filedata = new char[ file_size ];
         while ( thefData )
         { // Loop checking, 11.05.2015, T. Koi
            thefData.read( filedata , file_size );
         }
         thefData.close();
         data = new G4String ( filedata , file_size );
         delete [] filedata;
      }
      else
      {
         // found no data file
         // set error bit to the stream   
         iss.setstate( std::ios::badbit ); 
      }
   }
   if ( data != 0 )
   {
      iss.str(*data);
      G4String id;
      iss >> id;
      if ( id == "G4NDL" )
      {
         //Register information of file
         G4String source;
         iss >> source;
         register_data_file(filename,source);
      }
      else
      {
         iss.seekg( 0 , std::ios::beg );
      }
   }
   in->close(); delete in;
   delete data;
}

void G4ParticleHPManager::GetDataStream2( G4String filename , std::istringstream& iss ) 
{
   // Checking existance of data file 

   G4String compfilename(filename);
   compfilename += ".z";
   std::ifstream* in = new std::ifstream ( compfilename , std::ios::binary | std::ios::ate );
   if ( in->good() )
   {
      // Compressed file is exist 
      in->close(); 
   }
   else
   {
      std::ifstream thefData( filename , std::ios::in | std::ios::ate );
      if ( thefData.good() )
      {
         // Regular text file is exist
         thefData.close();
      }
      else
      {
         // found no data file
         // set error bit to the stream   
         iss.setstate( std::ios::badbit ); 
      }
   }
   delete in;
}

void G4ParticleHPManager::SetVerboseLevel( G4int newValue )
{
   G4cout << "You are setting a new verbose level for Particle HP package." << G4endl;
   G4cout << "the new value will be used in whole of the Particle HP package, i.e., models and cross sections for Capture, Elastic, Fission and Inelastic interaction." << G4endl;
   verboseLevel = newValue;
}

void G4ParticleHPManager::register_data_file(G4String filename, G4String source)
{
   mDataEvaluation.insert( std::pair < G4String , G4String > ( filename , source ) );
}

void G4ParticleHPManager::DumpDataSource()
{

   G4cout << "Data source of this Partile HP calculation are " << G4endl;
   for (  std::map< G4String , G4String >::iterator 
          it = mDataEvaluation.begin() ; it != mDataEvaluation.end() ; it++ )
   {
      G4cout << it->first << " " << it->second << G4endl;
   }
   G4cout << G4endl;
}

G4PhysicsTable* G4ParticleHPManager::GetInelasticCrossSections(const G4ParticleDefinition* particle )
{
   if ( theInelasticCrossSections.end() !=  theInelasticCrossSections.find( particle ) )
      return theInelasticCrossSections.find( particle )->second; 
   else 
      return 0; 
}

void G4ParticleHPManager::RegisterInelasticCrossSections( const G4ParticleDefinition* particle, G4PhysicsTable* val )
{
   theInelasticCrossSections.insert( std::pair<const G4ParticleDefinition* , G4PhysicsTable* >( particle , val ) ); 
}

std::vector<G4ParticleHPChannelList*>* G4ParticleHPManager::GetInelasticFinalStates(const G4ParticleDefinition* particle)
{
   if ( theInelasticFSs.end() != theInelasticFSs.find( particle ) )
      return theInelasticFSs.find( particle )->second;
   else 
      return 0;
}

void G4ParticleHPManager::RegisterInelasticFinalStates( const G4ParticleDefinition* particle , std::vector<G4ParticleHPChannelList*>* val )
{
   theInelasticFSs.insert ( std::pair<const G4ParticleDefinition*,std::vector<G4ParticleHPChannelList*>*>( particle , val ) ); 
}
