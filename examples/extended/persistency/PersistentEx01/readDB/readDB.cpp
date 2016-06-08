#include "HepODBMS/clustering/HepDbApplication.h"

#include "G4PEvent.hh"
#include "G4PPrimaryVertex.hh"
#include "G4PPrimaryParticle.hh"
#include "G4PHCofThisEvent.hh"
#include "G4PVHitsCollection.hh"

class dbAccessApp : public HepDbApplication {
  // Application inherits session control from HepDbApplication
public:
  // this application implements just one method: run the 
  dbAccessApp(const char *name) : HepDbApplication(name)
  {};
  
  int run()
  {
    // print an 
    message("about to initialise the db connection");
    Init();        // initialise the db session
    message("starting a read transaction");
    startRead(); // start a read transaction

    HepDatabaseRef  myDb = db("Events");

    // if the database ref is not valid:
    // - print a message
    // - exit the application with an error code 
    if (myDb == 0)
      fatal("could not find Events DB");
    
    // locates Events container in this database
    HepContainerRef cont = container("EventContainer"); 
    if (cont == 0 )
      fatal("could not find Events Container");

    // initialize iterator for G4PEvent in "Events" Container
    ooItr(G4PEvent) pevent_iterator;
    pevent_iterator.scan(cont);

    // Loop for all G4PEvent in the "EventContainer"

    while (pevent_iterator.next())
    {
      // access this G4PEvent
      int evt_id = pevent_iterator->GetEventID();
      int n_pvertex = pevent_iterator->GetNumberOfPrimaryVertex();

      cout << endl << "Reading event #" << evt_id << ":" << endl;
      cout << "  No. of primary vertex: " << n_pvertex << endl;

      // Loop for all primary vertex in this event
      for ( int i = 0; i < n_pvertex; i++ )
      {
        HepRef(G4PPrimaryVertex) pvertex = 
                             pevent_iterator->GetPrimaryVertex(i);
        cout << "    ObjectID of the primary vertex: "
             << pvertex.sprint() << endl;
        cout << "    No. of particle in the primary vertex: "
             << pvertex->GetNumberOfParticle() << endl;

        for ( int j = 0; j < pvertex->GetNumberOfParticle() ; j++ )
        {
           HepRef(G4PPrimaryParticle) particle = pvertex->GetPrimary(j);
           cout << "      PDGcode: " << particle->GetPDGcode()
                << "      Px = " << particle->GetPx()
                << "      Py = " << particle->GetPy()
                << "      Pz = " << particle->GetPz() << endl;
        }
      }

      // HitCollections of this event
      const HepRef(G4PHCofThisEvent) hcte = pevent_iterator->GetHCofThisEvent();
      if( hcte != 0 )
      {
        cout << "  ObjectID of the Hit Collections: " << hcte.sprint() << endl;

        // Capacity of the Hit Collections VArray
      //  G4int numHC = hcte->GetCapacity();

        // Actual number of Hit Collections in this event
        G4int numHCact = hcte->GetNumberOfCollections();

        cout << "  No. of the Hit Collections: " << numHCact << endl;
        for( G4int i = 0; i<numHCact; i++)
        {
          const HepRef(G4PVHitsCollection) aHC = hcte->GetHC(i);
          if( aHC != 0 )
          {
            cout << "    Hit Collection[" << i << "]:  Name: "
                 << aHC->GetName();
            cout << "     Sensitive Detector: " << aHC->GetSDname() << endl;
          }
        }
      }
      else
      { cout << "  No Hit Collections for this event." << endl; }
    }

    message("End of loop over events.");

    // finish this read transaction and return
    commit();
    return 0;
  }
  
};


int main(int argc, const char *argv[])
{
  dbAccessApp myApp(argv[0]);  // create an application object
  
  return myApp.run();    // call it's run method
}


