#ifndef G4TBBTASK_HH
#define G4TBBTASK_HH

#include <tbb/task.h>
#include "G4Types.hh"
#include "G4String.hh"

class G4VtbbJob;

class G4tbbTask : public tbb::task {
public:
  G4tbbTask(G4VtbbJob* job, 
            G4int i_event,  /* G4int n_select, const G4String& msg,*/ 
            G4int seedslenght);
  virtual ~G4tbbTask();
  tbb::task* execute();
private:
  //const G4int select;
  const G4int event;
  //const G4String msg;
  long* seeds;
  G4VtbbJob* job;
}; 
#endif //G4TBBTASK_HH
