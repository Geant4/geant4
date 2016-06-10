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
#include "G4MPIutils.hh"
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <assert.h>
#include <utility>
#include "globals.hh"



G4mpi::commMap_t G4mpi::buildCommunicationMap(
    std::vector<G4mpi::rank_t>& input ) {
  using namespace G4mpi;
  //Check validity of input
  std::sort(input.begin(),input.end());
  if ( input.size() < 1 || input[0] != 0 ) {
      G4Exception("G4mpi::buildCommunicationMap(...)","G4mpi001",FatalException,
          "Empty input or cannot find rank 0 in input.");
  }
  //Requested that no duplicates!
  std::vector<rank_t> copy(input.size());
  std::copy(input.begin(),input.end(),copy.begin());
  copy.erase( std::unique(copy.begin(),copy.end()),copy.end());
  if ( copy != input )
    {
      G4Exception("G4mpi::buildCommunicationMap(...)","G4mpi001",FatalException,
          "There are duplicates in list of input ranks.");
    }
  
  //The final communication map
  commMap_t mymap;
  //The communication map key
  int cycle = 0;
  //The communication map value
  std::vector<couple_t> couples;
  //Start a loop (on cycles) that will break
  do {
    //An helper container
    std::vector<rank_t> receiving;
    couples.clear();
    //Loop on all input until there is
    //at least a couple
    while ( input.size() > 1 ) {
      //Sort input in ascending order
      std::sort( input.begin(),input.end() );
      //Pop from back of the input a couple
      const auto& send_ = input.back();
      input.pop_back();
      const auto& rec_ = input.back();
      input.pop_back();
      //Actually the receiving is not empty,
      //remember it because we have to add it back
      receiving.push_back( rec_ );
      //This is a couple for this cycle
      couples.push_back( std::make_pair(send_,rec_) );
    }
    //Populate final map for this cycle
    mymap[cycle++]=couples;
    //Let's put back in the input container the receivers
    input.insert( input.end() , receiving.begin() , receiving.end() );
    //Let's continue until there is ony one rank in input (number 0)
  } while ( input.size()!=1 );
  return mymap;
}

//Test function
int _testMe(int argc,char** argv) {
  using namespace G4mpi;
  unsigned int worldSize = 10;
  if ( argc > 1 ) worldSize = atoi(argv[1]);
  unsigned int myRank = worldSize-1;
  if ( argc > 2 ) myRank = atoi(argv[2]);
  std::cout<<"World size: "<<worldSize<<std::endl;

  assert( myRank < worldSize);
  //MPI function stubs
  auto MPI_Receive = [](const rank_t& s, const rank_t& r) {
    std::cout<<"MPI_Receive from: "<<s<<" to "<<r<<std::endl;
    return 0;
  };
  auto MPI_Send = [](const rank_t& s, const rank_t& r) {
    std::cout<<"MPI_Send from: "<<s<<" to "<<r<<std::endl;;
    return 0;
  };
  auto MPI_Barrier = [] {
    std::cout<<"MPI_Barrier"<<std::endl;
    return 0;
  };

  //Build the initial network of ranks
  std::vector<rank_t> ranks(worldSize);
  for ( unsigned int i = 0 ; i<worldSize ; ++i ) {
    if ( i != 2 )
      { ranks.push_back(i);
      } else {
      ranks.push_back(i);
      ranks.push_back(i);
    }
  }
  //Remove duplicates
  ranks.erase( std::unique(ranks.begin(),ranks.end()),ranks.end());
  
  //Optimize network trafic
  if ( ranks.size() == 1 ) {
    std::cout<<"only one rank, nothing to do"<<std::endl;
    return 0;
  }
  auto comms = G4mpi::buildCommunicationMap( ranks );
  assert( ranks.size() == 1 && ranks[0] == 0 );

  std::cout<<"Communiction Map (size: "<<comms.size()<<"):"<<std::endl;
  for ( const auto& x : comms ) {
    std::cout<<"Cycle "<<x.first<<": ";
    for ( const auto& y : x.second ) {
      std::cout<<y.first<<"->"<<y.second<<", ";
    }
    std::cout<<std::endl;
  }

  std::cout<<"Simulate communication pattern for rank: "<<myRank<<std::endl;
  for (const auto& x: comms ) {
    std::cout<<"Cycle "<<x.first<<std::endl;
    for ( const auto& y : x.second ) {
      if ( myRank == y.first ) { MPI_Send(y.first,y.second); }
      else if ( myRank == y.second ) { MPI_Receive(y.first,y.second); }
    }
    //Important: Wait for this cycle to end before going to the next, even if
    //this rank did not do anything
    //This is needed to be sure that the redcutions are done correctly
    MPI_Barrier();
  }
  return 0;
}

void G4mpi::Merge( std::function<void(unsigned int)> senderF ,
                                    std::function<void(unsigned int)> receiverF,
                                    std::function<void(void)> barrierF ,
                                    unsigned int commSize ,
                                    unsigned int myrank) {
      //Optimize communications between ranks
      std::vector<G4mpi::rank_t> ranks(commSize);
      std::iota(ranks.begin(),ranks.end(),0); //{0,1,2,3,...}
      auto comms = G4mpi::buildCommunicationMap(ranks);
      //Loop on cycles of communications
      for ( const auto& cycle : comms ) {
          //Each cycle is a set of communications between ranks, it is guarantted that
          //each rank participate in one and only one communication for each cycle
         for (const auto& pattern : cycle.second ) {
             //pattern is a couple: sender,receiver
             if ( myrank == pattern.first ) {
                 //Send to destination
                 senderF(pattern.second);
             }
             else if ( myrank == pattern.second ) {
                 //Receive from source
                 receiverF(pattern.first);
             }
         }
         //Important: Wait for this cycle to end before going to the next, even if this rank
         //did not do anything
          //This is needed to be sure that the redcutions are done correctly
         barrierF();
    }
}
