//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//The template for split classes as follows:
//G4LogicalVolume, G4Region, G4VPhysicalVolume, G4PolyconeSide
//G4PolyhedraSide, G4PVReplica, G4Material. 
#ifndef G4MTPRIVATESUBINSTANCEMANAGERER_HH
#define G4MTPRIVATESUBINSTANCEMANAGERER_HH

#include <stdlib.h>
#include <string.h>

template <class G4MTPrivateObject>
class G4MTPrivateSubInstanceManager
{
public:
  static __thread G4MTPrivateObject* offset;

  G4MTPrivateSubInstanceManager() {
    totalobj = 0;
    totalspace = 0;
    sharedOffset = 0;
  };

  //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  //Invoked by the master thread to create a new subinstance whenever
  //a new split class instance is created.
  int CreateSubInstance() {
    totalobj++;
    if (totalobj > totalspace)
    {
      totalspace=totalspace + 512; //DYNAMICINCREASESTEP                                                     
      offset = (G4MTPrivateObject *) realloc(offset, totalspace * sizeof(G4MTPrivateObject));
      if (offset == NULL)
      {
        printf("Can not malloc space in G4MTPrivateSubInstanceManager\n");
        exit(-1);
      }
      sharedOffset = offset;
    }
    return (totalobj - 1);
  };

  //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  //Invoked by each worker thread to copy all the subinstance array
  //from the master thread.
  void SlaveCopySubInstanceArray()
  {
    if (offset) return;
    offset = (G4MTPrivateObject *) realloc(offset, totalspace * sizeof(G4MTPrivateObject));
    if (offset == NULL)
    {
      printf("Can not malloc space in G4MTPrivateSubInstanceManager\n");
      exit(-1);
    }
    memcpy(offset, sharedOffset, totalspace * sizeof(G4MTPrivateObject));
  }

  //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  //Invoked by each worker thread to create the subinstance array and
  //initialize each subinstance using a particular method defined by
  //the subclass.
  void SlaveInitializeSubInstance()
  {
    if (offset) return;
    offset = (G4MTPrivateObject *) realloc(offset, totalspace * sizeof(G4MTPrivateObject));

    if (offset == NULL)
    {
      printf("Can not malloc space in G4MTPrivateSubInstanceManager\n");
      exit(-1);
    }

    for (int i = 0 ; i < totalspace ; i++)
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
  int totalobj;
  int totalspace;
  G4MTPrivateObject* sharedOffset;
};

#endif
