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
#include "G4MPIScorerMerger.hh"
#include <map>
#include <strstream>

G4MPIScorerMerger::G4MPIScorerMerger( G4ScoringManager* mgr, 
                                      G4int destination,
                                      G4int verbosity ) : 
scoringManager(mgr),commSize(0),destinationRank(destination),verbose(verbosity)
{}

#define DMSG( LVL , MSG ) { if ( verbose > LVL ) { G4cout << MSG << G4endl; } }

std::ostream& operator<<(std::ostream& os , 
                         const G4MPIScorerMerger::convMap_t& cnv ) {
  static const G4int maxelems = 10;
  os<<" Name: "<<cnv.name<<" with : "<<cnv.numElems<<" elements\n";
  os<<"\tIndexes :";
  for ( G4int i = 0 ;
        i < ((cnv.numElems<maxelems) ? cnv.numElems : maxelems) ;
          ++i ) os<<" "<<(cnv.indexes)[i];
  if ( cnv.numElems>maxelems ) os<<" ...";
  for ( G4int i = ( (cnv.numElems-maxelems > 0) ? cnv.numElems-maxelems : cnv.numElems );
        i < cnv.numElems ; ++i ) os<<" "<<(cnv.indexes)[i];
  os<<"\n\tValues :";
  for ( G4int i = 0 ;
        i < ((cnv.numElems<maxelems) ? cnv.numElems : maxelems) ;
          ++i ) os<<" "<<(cnv.values)[i];
  if ( cnv.numElems>maxelems ) os<<" ...";
  for ( G4int i = ( (cnv.numElems-maxelems > 0) ? cnv.numElems-maxelems : cnv.numElems );
        i < cnv.numElems ; ++i ) os<<" "<<(cnv.values)[i];
  return os;
}

std::ostream& operator<<(std::ostream& os, 
                         G4THitsMap<double>& map) 
{
  os<<map.GetName()<<" "<<map.GetSDname()<<" "<<map.GetSize()<<"\n";
  for ( std::map<G4int,G4double*>::const_iterator it =
        map.GetMap()->begin(); it != map.GetMap()->end() ;
        ++it)
    os<<it->first<<" "<<*(it->second)<<"\n";
  return os;
}

G4MPIScorerMerger::convMap_t*
G4MPIScorerMerger::convertMap( const G4String& mapName ,
                               G4THitsMap<double>* map ) const 
{
  DMSG( 2 , "Converting G4THitsMap<double> "<<map<<
       " with name "<<mapName);
  convMap_t* converted = new convMap_t;
  converted->name = mapName;
  DMSG(2,converted->name);
  converted->numElems = map->GetSize();
  DMSG(2,converted->numElems);
  converted->indexes = new G4int[converted->numElems];
  converted->values = new G4double[converted->numElems];
  std::map<G4int,double*>* mm = map->GetMap();
  G4int counter=0;
  for ( std::map<G4int,G4double*>::const_iterator it = mm->begin();
        it != mm->end() ; ++it ) {
    //DMSG(2,it->first<<" "<<*(it->second)<<" "<<counter);
    (converted->indexes)[counter] = it->first;
    (converted->values)[counter++] = *(it->second);
    }
  DMSG( 2 , "Converted to: "<<*converted );
  return converted;
}

void
G4MPIScorerMerger::convertMesh( const G4VScoringMesh* mesh )
{
  DMSG(2,"Coverting G4VScoringMesh: "<<mesh);
  clear();
  const MeshScoreMap& map = mesh->GetScoreMap();
  DMSG(2,"Converting "<<map.size()<<" score maps");
  for ( MeshScoreMap::const_iterator it = map.begin() ;
        it != map.end() ; ++it ) 
    {
      convertedMesh.push_back( convertMap(it->first,it->second) );
    }
  DMSG(2,"Conversion of G4VScoringMesh: "<<mesh<<" done");
}

void 
G4MPIScorerMerger::clear() 
{
  for ( std::vector<convMap_t*>::iterator it = convertedMesh.begin() ;
        it != convertedMesh.end() ; ++it )
    {
      delete[] (*it)->indexes;
      delete[] (*it)->values;
      delete *it;
      *it = 0;
    }
  convertedMesh.erase(convertedMesh.begin(),convertedMesh.end());
}


void G4MPIScorerMerger::Merge()
{
  DMSG(0,"G4MPIScorerMerger::Merge() called");
  G4int myrank =  MPI::COMM_WORLD.Get_rank();
  commSize = MPI::COMM_WORLD.Get_size();
  COMM_G4COMMAND_ = MPI::COMM_WORLD.Dup();
  DMSG(0,"Comm world size: "<<commSize<<" this rank is: "
       <<myrank<<" sending to rank "<<destinationRank
       <<" Number of mesh: "<< scoringManager->GetNumberOfMesh() );
  for ( size_t i = 0 ; i < scoringManager->GetNumberOfMesh() ; ++i )
    {
      if ( myrank != destinationRank ) {
        meshID = static_cast<G4int>(i);
        SendOneMesh();
      } else { 
        ReceiveOneMesh();
      }
    }
  DMSG(0,"G4MPIScorerMerger::Merge done");
}

void G4MPIScorerMerger::SendOneMesh() 
{
  DMSG(1,"Sending mesh with ID: "<<meshID);
  G4VScoringMesh* mesh = scoringManager->GetMesh(meshID);
  convertMesh( mesh );
  COMM_G4COMMAND_.Send(&meshID,1,MPI::INT,
                       destinationRank,G4MPImanager::kTAG_DATA);
  G4int numelems=convertedMesh.size();
  DMSG(2,"Sending "<<numelems<<" maps");
  COMM_G4COMMAND_.Send(&numelems,1,MPI::INT,
                       destinationRank,G4MPImanager::kTAG_DATA);
  for ( std::vector<convMap_t*>::const_iterator it = convertedMesh.begin() ;
        it != convertedMesh.end() ; ++it ) {
    const convMap_t* elem = *it;
    DMSG(2,"Sending map: "<<elem);
    COMM_G4COMMAND_.Send(elem->name.c_str(),elem->name.length(),MPI::CHAR,
                         destinationRank,G4MPImanager::kTAG_DATA);
    COMM_G4COMMAND_.Send(&(elem->numElems),1,MPI::INT,
                         destinationRank,G4MPImanager::kTAG_DATA);
    COMM_G4COMMAND_.Send(elem->indexes,elem->numElems,MPI::INT,
                         destinationRank,G4MPImanager::kTAG_DATA);
    COMM_G4COMMAND_.Send(elem->values,elem->numElems,MPI::DOUBLE,
                         destinationRank,G4MPImanager::kTAG_DATA);
  }
  DMSG(1,"Sending of mesh with ID: "<<meshID<<" Done.");
}

void G4MPIScorerMerger::ReceiveOneMesh() 
{
  DMSG(1,"Receiving of mesh");
  clear();
  for ( G4int rank = 0 ; rank<commSize; ++rank ) {
    if ( rank == destinationRank ) continue; // Do not receive from myself
    COMM_G4COMMAND_.Recv(&meshID,1,MPI::INT,rank,G4MPImanager::kTAG_DATA);
    G4int numElems = 0;
    COMM_G4COMMAND_.Recv(&numElems,1,MPI::INT,rank,G4MPImanager::kTAG_DATA);
    for ( G4int i = 0 ; i<numElems ; ++i ) {
      convMap_t* elem = new convMap_t;
      //G4int strlen = 0;
      MPI::Status status;
      COMM_G4COMMAND_.Probe(rank,G4MPImanager::kTAG_DATA,status);
      G4int strlen = status.Get_count(MPI::CHAR);
      char* buf = new char[strlen];
      COMM_G4COMMAND_.Recv(buf,strlen,MPI::CHAR,rank,
                           G4MPImanager::kTAG_DATA,status);
      elem->name = G4String(buf,strlen);
      delete[] buf;
      COMM_G4COMMAND_.Recv(&(elem->numElems),1,MPI::INT,rank,
                           G4MPImanager::kTAG_DATA);
      elem->indexes = new G4int[elem->numElems];
      elem->values = new G4double[elem->numElems];
      COMM_G4COMMAND_.Recv(elem->indexes,elem->numElems,MPI::INT,rank,
                           G4MPImanager::kTAG_DATA);
      COMM_G4COMMAND_.Recv(elem->values,elem->numElems,MPI::DOUBLE,rank,
                           G4MPImanager::kTAG_DATA);
      convertedMesh.push_back(elem);
      DMSG(2,"Received one mesh map: "<<*elem);
    }
    DMSG(2,"Received one mesh, with "<<convertedMesh.size()<<" maps");
    //Received from Rank number rank
    MergeOneMesh();
  }
  DMSG(1,"Receiving of mesh done");
}

void G4MPIScorerMerger::MergeOneMesh()
{
  DMSG(2,"Merging one mesh");
  G4VScoringMesh* mesh = scoringManager->GetMesh(meshID);
  if ( ! mesh ) {
    G4ExceptionDescription msg;
    msg<<"Cannot find mesh with id: "<<meshID;
    G4Exception("G4MPIScorerMerger::MergeOneMesh()","G4MPI001",FatalException,
                 msg);
  }
  for ( std::vector<convMap_t*>::const_iterator it = convertedMesh.begin() ; 
        it != convertedMesh.end() ; ++it ) 
  {
    //Create a hits-collection from this convMap_t object
    const convMap_t* const elem = *it;
    G4THitsMap<G4double> hc(mesh->GetWorldName(),elem->name);
    for ( G4int i = 0 ; i < elem->numElems; ++i )
      hc.set( elem->indexes[i] , elem->values[i] );
    DMSG(3,"Original mesh: "<<*(mesh->GetScoreMap().find(elem->name)->second));
    DMSG(3,"Received mesh: "<<hc);
    mesh->Accumulate(&hc);
    DMSG(3,"Original mesh after accumulation: "
         <<*(mesh->GetScoreMap().find(elem->name)->second));
  }
  DMSG(2,"Merging one mesh done");
}

