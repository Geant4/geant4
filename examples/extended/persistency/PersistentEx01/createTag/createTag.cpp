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

    message("starting an update transaction");
    startUpdate(); // start an update transaction

    // create an iterator over all events in the container "Events" in db "Events"
    db("Events");
    ooItr(G4PEvent) pevent_iterator;
    pevent_iterator.scan(container("EventContainer"));

#ifdef REMOVE_OLD_TAGS
    // remove an existing tag collection with the name <tagName>
    HepExplorable *cd = HepExplorable::findExplorable(tagName);
    if (cd) {
      cout << "removing old tag:" << tagName << endl;
      cd->removeDescription();
      delete cd;
    }
#endif    

    HepExplorableGenericTags events; // create a tag collection ...
    db("Tags");                      // ... in db "Tags"

    if (events.createDescription(tagName) == 0) // start creating a new tag field description
      fatal("could not create new tag");  // print the message and stop the application

 
    // define some attributes in the tag
    TagAttribute<long>   event_id(events,"eventID");
    TagAttribute<long>   n_pvertex(events,"NumberOfPrimaryVertex");
    
    
    // Loop over all tag collection elements (= all events)

    while(pevent_iterator.next())
    {
      HepRef(G4PEvent) evt = pevent_iterator;
      // create a new tag entry and bind it to the current event
      events.newTag(evt);
      
      // fill the tag variables with data from the event
      event_id = evt->GetEventID();
      n_pvertex = evt->GetNumberOfPrimaryVertex();
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

