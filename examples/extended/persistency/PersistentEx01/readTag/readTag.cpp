#include "HepODBMS/tagdb/HepTagDbApplication.h"
#include "HepODBMS/tagdb/HepAnalysis.h"

#include "G4PEvent.hh"

const char tagName[] = "geant4Events";

class dbAccessApp : public HepTagDbApplication {
  // Application inherits session control from HepTagDbApplication
public:
  // this application implements just one method: run the 
  dbAccessApp(const char *name) : HepTagDbApplication(name)
  {};
  
  int run()
  {
    // print an 
    message("about to initialise the db connection");
    Init();        // initialise the db session
    message("starting a read transaction");
    startRead(); // start a read transaction

    // find the tag collection by name
    HepExplorable *events =
      HepExplorable::findExplorable(tagName);
    
    if (events == 0)
       fatal("could not find tag collection");
 
    // define some attributes in the tag
    TagAttribute<long>   event_id(events,"eventID");
    TagAttribute<long>   n_pvertex(events,"NumberOfPrimaryVertex");
    
    // also use G4PEvent directly...
    HepRef(G4PEvent) g4evt;
    
    // Loop over all tag collection elements (= all events)

    for(events->start(); events->next() != 0; )
    {
      // use a tag attribute directly in a print statment
      cout << "Reading event #" << event_id << ":" << G4endl;
      cout << "  number of primary vertex: " << n_pvertex << G4endl;
      
      //
      // [ could do some selection based on tag attributes here ... ]
      //
      //  eg.   if( n_pvertex < 20 )
      //        { ...do something... }
      //        else
      //        { ...do anything else... }
      //

      // retrieve the associated event object for the current tag
      // and navigate to the associated objects/members which are
      // not present in the tag (GetNumberOfParticle() in this example).

      if (getEvent(g4evt, events) )
      {
        cout << "  ObjectID of event: " << g4evt.sprint() << G4endl;
        cout << "  No. of the primary vertex (from G4PEvent): "
             << g4evt->GetNumberOfPrimaryVertex() << G4endl;
  
        for ( int i = 0; i < n_pvertex; i++ )
        {
          HepRef(G4PPrimaryVertex) pvertex = g4evt->GetPrimaryVertex(i);
          cout << "    ObjectID of primary vertex: " << pvertex.sprint() << G4endl;
          cout << "    No. of particle in the primary vertex: "
               << pvertex->GetNumberOfParticle() << G4endl;
        }
      }
      else
      { cout << "  No event object found." << G4endl; }

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

