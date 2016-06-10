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
// Utility functions for MPI G4 interface
#ifndef G4MPIUTILS_HH
#define G4MPIUTILS_HH
#include <map>
#include <vector>
#include <functional>
#include <numeric>

//Namespace with some utility functions for G4 MPI integration.
//Main utilities:
// G4mpi::Merge(...) : Merge results via a semi-optimized communication
//                                patterns between ranks. Note that this implementation
//                                is not topology aware. This means that MPI_Reduce and
//                                MPI_Gather are more performant. However if you cannot
//                                implement an appropriate MPI reducer or you cannot efford
//                                the memory overhead of Gather, this can be used instead of
//                                p2p communications.
namespace G4mpi {
  //Simple data type representing a rank
  typedef unsigned int rank_t;
  //A couple of sending/receiving ranks
  typedef std::pair<rank_t,rank_t> couple_t;
  //This map represent, for each cycle (key) a set of communications
  //pairs
  typedef std::map<int,std::vector<couple_t> > commMap_t;

  //This function takes as input a vector of rank_t objects representing
  //a communication node identified by an ID (rank:int).
  //It returns a map of cycle:int -> vector<pairs<rank_t> >
  //Representing a sequence of communciation cycles. At each communication cycle
  //one or more p2p communications are established: in the pair the first element
  //is the sender and the second element of the pair is the receiver
  //At the end of the cycles all communications have been done to rank 0
  //For example with 4 nodes: [0,1,2,3] we have:
  //   Cycle 0: (3->2),(1->0)
  //   Cycle 1: (2->0)
  // With 5 nodes:
  //   Cycle 0: (4->3),(2->1)
  //   Cycle 1: (3->1)
  //   Cycle 2: (1->0)
  // The algorithm can be used to implement a communication across mpi ranks
  // optimizing the network trafic. Each rank (the nodes) can send/receive to another node.
  // Once they have sent out the payload they become empty and non-active anymore.
  commMap_t buildCommunicationMap( std::vector<rank_t>& input );

  //Performs merging to rank 0 using the provided sender, receiver and barrier functions.
  //CommSize is the size of the communicator and myrank is the rank of the caller
  //For example: assume a class UserMerger has two members Send(uint) and
  //                      Receive(uint) and we are using a MPI::Intracomm object as
  //                      communicator, then to use this function the ranks can:
  //                            using std::placeholers::_1;
  //                            std::function<void(unsigned int)> sender =
  //                                            std::bind(&Merger::Send,&mergerInst,_1);
  //                            std::function<void(unsigned int)> receiver =
  //                                            std::bind(&Merger::Receiver,&mergerInst,_1);
  //                            std::function<void(void)> barrier =
  //                                            std::bind(&MPI::Intracomm::Barrier,&commInst);
  //                            G4mpi::Merge(sender,receiver,barrier,commSize,myrank);
  void Merge( std::function<void(unsigned int)> senderF ,
                          std::function<void(unsigned int)> receiverF ,
                          std::function<void(void)> barrierF ,
                          unsigned int commSize , unsigned int myrank);
  //Type representing a merging functions
  typedef std::function<void(std::function<void(unsigned int)>,
                             std::function<void(unsigned int)>,
                             std::function<void(void)>,
                             unsigned int, unsigned int)>
          mergerHandler_t;
}

#endif //G4MPIUTILS_HH
