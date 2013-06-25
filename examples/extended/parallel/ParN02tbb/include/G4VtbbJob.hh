#ifndef G4VTBBJOB_HH
#define G4VTBBJOB_HH

#include "G4String.hh"
class G4tbbRunManager;

class G4VtbbJob {
  friend class G4tbbTask;
public:
  G4VtbbJob(const G4String& macro="");
  virtual ~G4VtbbJob();
  virtual void InitRun( G4tbbRunManager* rm );
protected:
  virtual void UserActions(G4tbbRunManager* rm=0 ) =0;
  virtual void InitSetup(G4tbbRunManager* rm=0 ) =0;
  virtual void JobPrepare(G4tbbRunManager* rm=0 ) =0;
private:
  void ThreadSafeInitSetup(G4tbbRunManager* rm);
  G4String macroFile;
};

#endif //G4VTBBJOB_HH
