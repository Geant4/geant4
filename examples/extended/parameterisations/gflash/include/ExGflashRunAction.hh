#ifndef ExGflashRunAction_h
#define ExGflashRunAction_h
using namespace std;
#include "G4UserRunAction.hh"
#include "globals.hh"


class G4Run;

class ExGflashRunAction : public G4UserRunAction
{
  public:
  
  
  ExGflashRunAction();
  ~ExGflashRunAction();
  void BeginOfRunAction(const G4Run* aRun);
  void EndOfRunAction(const G4Run* aRun);
  void setRunID (int i) {runID=i;};
  private:
    	G4int runID;
};


#endif
