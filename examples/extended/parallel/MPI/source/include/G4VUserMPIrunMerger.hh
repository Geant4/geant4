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
#ifndef G4MPIRUNMERGER_HH
#define G4MPIRUNMERGER_HH
#include "G4Run.hh"
#include <mpi.h>
#include "G4MPImanager.hh"

class G4VUserMPIrunMerger {
public:
  G4VUserMPIrunMerger();
  G4VUserMPIrunMerger( const G4Run* aRun ,
                  G4int destination = G4MPImanager::kRANK_MASTER ,
                  G4int verbosity = 0);
  virtual ~G4VUserMPIrunMerger() { if ( ownsBuffer) DestroyBuffer(); }
  void SetRun( G4Run* r ) { run = r; }
  void SetDestinationRank( G4int i ) { destinationRank = i; }
  void SetVerbosity( G4int ver ) { verbose = ver; }

  virtual void Merge();

protected:
  virtual void Pack() = 0;
  virtual G4Run* UnPack() = 0;

  void InputUserData( /*const*/ void* input_data ,const  MPI::Datatype& dt, int count) {
    input_userdata.push_back( const_registered_data{input_data,dt,count} );
  }
  void OutputUserData( void* input_data ,const  MPI::Datatype& dt, int count) {
    output_userdata.push_back( registered_data{input_data,dt,count} );
  }

 // void GetUserData(void* output_data,const MPI::Datatype& dt, int count);

  void SetupOutputBuffer(char* buff, G4int size, G4int position) {
    outputBuffer = buff;
    outputBufferSize=size;
    outputBufferPosition=position;
  }
  void DestroyBuffer() {
    delete[] outputBuffer;
    outputBuffer = nullptr;
    outputBufferSize=0;
    outputBufferPosition=0;
    ownsBuffer = false;
  }

  G4int GetPosition() const { return outputBufferPosition; }
  char* GetBuffer() const { return outputBuffer; }
  G4int GetBufferSize() const { return outputBufferSize; }

  void Send(const unsigned int destination);
  void Receive(const unsigned int source);
private:
  char* outputBuffer;
  G4int outputBufferSize;
  G4int outputBufferPosition;
  G4bool ownsBuffer;
  unsigned int destinationRank;
  G4Run* run;
  unsigned int commSize;
  MPI::Intracomm COMM_G4COMMAND_;
  G4int verbose;
  long bytesSent;

  //Input data to send (read-only)
  struct const_registered_data {
    const_registered_data(const const_registered_data&) = default;
    const_registered_data& operator=(const const_registered_data&) = default;
    //const_registered_data(const_registered_data&&) = default;
    //const_registered_data& operator=(const_registered_data&&) = default;
    /*const*/ void* p_data;
    /*const*/ MPI::Datatype dt;
    /*const*/ int count;
  };
  std::vector<const_registered_data> input_userdata;

  //Output data
    struct registered_data {
    registered_data(const registered_data&) = default;
    registered_data& operator=(const registered_data&) = default;
    void* p_data;
    /*const*/ MPI::Datatype dt;
    /*const*/ int count;
  };
    std::vector<registered_data> output_userdata;

};


#endif //G4MPIRUNMERGER_HH

