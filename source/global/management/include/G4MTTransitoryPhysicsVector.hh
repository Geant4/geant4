//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//The template for split classes as follows:
//G4PhysicsVector
#ifndef G4MTPRIVATEPHYSICSVECTOR_HH
#define G4MTPRIVATEPHYSICSVECTOR_HH

#include <stdlib.h>
#include <string.h>

#ifdef G4MULTITHREADED
extern pthread_mutex_t mutexPhysicsVector;
//pthread_mutex_t mutexPhysicsVector = PTHREAD_MUTEX_INITIALIZER;
#endif

template <class G4MTPrivateObject>
class G4MTPrivatePhysicsVectorCounter
{
public:
  static G4MTPrivateObject* offsetshadow;
  static G4ThreadLocal G4MTPrivateObject* offset;

  G4MTPrivatePhysicsVectorCounter() {
  };

  //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  //Invoked by the master or work thread to create a new subinstance
  //whenever a new split class instance is created.
  int CreateSubInstance() {
#ifdef G4MULTITHREADED
    pthread_mutex_lock(&mutexPhysicsVector);
#endif

    totalobj++;
    if (totalobj > totalspace) NewSubInstances();

    if (phaseshadow == 0)
    {
      totalobjshadow = totalobj;
      offsetshadow = offset;
    }

    int totalobjlocal = totalobj;

#ifdef G4MULTITHREADED
    pthread_mutex_unlock(&mutexPhysicsVector);
#endif

    return (totalobjlocal - 1);
  };

  //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  //Invoked by each worker thread to create the subinstance array and
  //initialize each subinstance using a particular method defined by
  //the subclass.
  void SlaveInitializeSubInstance()
  {
#ifdef G4MULTITHREADED
    pthread_mutex_lock(&mutexPhysicsVector);
#endif

    phaseshadow = 1;
    totalobj = totalobjshadow;
    totalspace = totalobj;
    offset = (G4MTPrivateObject *) malloc(totalspace * sizeof(G4MTPrivateObject));
    if (offset == NULL)
    {
      printf("Can not malloc space in G4MTPrivatePhysicsVectorCounter\n");
      exit(-1);
    }

    //    memcpy(offset, offsetshadow, totalspace * sizeof(G4MTPrivateObject));

    for (int i = 0 ; i < totalspace ; i++)
    {
      offset[i].initialize();
    }

#ifdef G4MULTITHREADED
    pthread_mutex_unlock(&mutexPhysicsVector);
#endif
  }

  //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  //Invoked by each worker thread to grow the subinstance array and
  //initialize each new subinstance using a particular method defined
  //by the subclass.
  void NewSubInstances()
  {
    if (totalspace  >= totalobj) return;
    int originaltotalspace = totalspace;
    totalspace = totalobj + 512; //DYNAMICINCREASESTEP
    offset = (G4MTPrivateObject *) realloc(offset, totalspace * sizeof(G4MTPrivateObject));
    if (offset == NULL)
    {
      printf("Can not malloc space in G4MTPrivatePhysicsVectorCounter\n");
      exit(-1);
    }

    for (int i = originaltotalspace ; i < totalspace ; i++)
    {
      offset[i].initialize();
    }
  }

  //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  //Invoked by all threads to free the subinstance array.
  void FreeSlave()
  {
    if (!offset) return;
    delete offset;
  }

private:
  static int phaseshadow; //0: master initialization, 1: worker initialization
  static int totalobjshadow;
  static G4ThreadLocal int totalobj;
  static G4ThreadLocal int totalspace;
};

#endif
