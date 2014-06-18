#include "G4MPIRunMerger.hh"
#include <mpi.h>

G4MPIRunMerger::G4MPIRunMerger( const G4Run* aRun , 
                                G4int destination ,
                                G4int ver) : 
   destinationRank(destination),
   run(const_cast<G4Run*>(aRun)),
   commSize(0),
   verbose(ver) {}

#define DMSG( LVL , MSG ) { if ( verbose > LVL ) { G4cout << MSG << G4endl; } }

void G4MPIRunMerger::SendDouble( G4double* val, G4int size ) 
{
  DMSG( 2 , "Sending double from "<<val<<" with size: "<<size);
  COMM_G4COMMAND_.Send(  val, size , MPI::DOUBLE, 
                       destinationRank,G4MPImanager::kTAG_DATA);
  DMSG( 2 , "Sent "<<( size > 1 ? val[0] : *val) );
}

void G4MPIRunMerger::SendInt( G4int* val , G4int size ) 
{
  DMSG( 2 , "Sending int from "<<val<<" with size: "<<size);
  COMM_G4COMMAND_.Send( val, size , MPI::INT, 
                       destinationRank,G4MPImanager::kTAG_DATA);
  DMSG( 2 , "Sent "<<( size > 1 ? val[0] : *val) );
}

void G4MPIRunMerger::ReceiveDouble( G4int rank, G4double* val , G4int size )
{
  DMSG( 2 , "Receiving double at "<<val<<" with size "<<size );
  COMM_G4COMMAND_.Recv( val, size, MPI::DOUBLE, rank , G4MPImanager::kTAG_DATA);
  DMSG( 2 , "Received "<<( size > 1 ? val[0] : *val) );
}

void G4MPIRunMerger::ReceiveInt( G4int rank, G4int* val , G4int size )
{
  DMSG( 2 , "Receiving int at "<<val<<" with size "<<size );
  COMM_G4COMMAND_.Recv( val, size, MPI::INT, rank , G4MPImanager::kTAG_DATA);
  DMSG( 2 , "Received "<<( size > 1 ? val[0] : *val) );
}

void G4MPIRunMerger::Send() 
{
  G4int nevts = run->GetNumberOfEvent();
  DMSG( 1 , "G4MPIRunMerger::Send() : Sending a G4run ("
        <<run<<") with "<<nevts<<" events,");
  SendInt( &nevts );
  DMSG( 1 , "G4MPIRunMerger::Send() : Done ");
}
void G4MPIRunMerger::Receive(G4int rank) 
{
  DMSG( 1 , "G4MPIRunMerger::Receive(...) : Receiving from rank "<<rank);
  G4Run* anEmptyRun = new G4Run;
  G4int nevts = 0;
  ReceiveInt( rank, &nevts );

  //Increment internal counter up to nevets
  for ( G4int i = 0 ; i<nevts ; ++i ) anEmptyRun->RecordEvent( NULL );

  //User data can go here
  //
  //Now merge received MPI run with global one
  DMSG(2,"Before G4Run::Merge : "<<run->GetNumberOfEvent());

  run->Merge( anEmptyRun );
  DMSG(2,"After G4Run::Merge : "<<run->GetNumberOfEvent());
  delete anEmptyRun;
}

void G4MPIRunMerger::Merge() 
{
  DMSG(0, "G4MPIRunMerger::Merge called");
  G4int myrank = MPI::COMM_WORLD.Get_rank();
  commSize = MPI::COMM_WORLD.Get_size();
  COMM_G4COMMAND_ = MPI::COMM_WORLD.Dup();
  DMSG(0,"Comm world size: "<<commSize<<" this rank is: "
       <<myrank<<" sending to rank "<<destinationRank);
  for ( G4int i = 0 ; i < commSize; ++i ) {
    //Send for all ranks except receiver
    if ( myrank != destinationRank ) Send();
  }
  //Receiver receives from all ranks
  if ( myrank == destinationRank ) {
    for (G4int fromRank = 0 ; fromRank < commSize; ++fromRank) {
      //Do not receive from myself
      if ( fromRank != destinationRank ) Receive(fromRank);
    }
  }
  DMSG(0,"G4MPIRunMerger::Merge done");
}

