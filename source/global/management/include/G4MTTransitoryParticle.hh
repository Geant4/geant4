//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//The template for split classes as follows:
//G4ParticleDefinition, G4VDecayChannel
#ifndef G4MTPRIVATEPARTICLECOUNTER_HH
#define G4MTPRIVATEPARTICLECOUNTER_HH

#include <stdlib.h>
#include <string.h>

template <class G4MTPrivateObject>
class G4MTPrivateParticleCounter
{
public:
  static G4ThreadLocal int slavetotalspace;
  static G4ThreadLocal G4MTPrivateObject* offset;

  G4MTPrivateParticleCounter() {
    totalobj = 0;
  };

  //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  //Invoked by the master or work thread to create a new subinstance
  //whenever a new split class instance is created. For each worker
  //thread, ions are created dynamically
  int CreateSubInstance() {
    totalobj++;
    if (totalobj > slavetotalspace) NewSubInstances();
    return (totalobj - 1);
  };

  //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  //Invoked by each worker thread to grow the subinstance array and
  //initialize each new subinstance using a particular method defined
  //by the subclass.
  void NewSubInstances()
  {
    if (slavetotalspace  >= totalobj) return;
    int originaltotalspace = slavetotalspace;
    slavetotalspace = totalobj + 512; //DYNAMICINCREASESTEP                                                                      
    offset = (G4MTPrivateObject *) realloc(offset, slavetotalspace * sizeof(G4MTPrivateObject));

    if (offset == NULL)
    {
      printf("Can not malloc space in G4MTPrivateParticleCounter\n");
      exit(-1);
    }

    for (int i = originaltotalspace ; i < slavetotalspace ; i++)
      offset[i].initialize();
  }

  //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  //Invoked by all threads to free the subinstance array.
  void FreeSlave()
  {
    if (!offset) return;
    delete offset;
  }

private:
  int totalobj;
};

#endif
