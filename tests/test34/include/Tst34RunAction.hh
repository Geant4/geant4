#ifndef Tst34RunAction_h
#define Tst34RunAction_h
using namespace std;
#include "G4UserRunAction.hh"
#include "globals.hh"


class G4Run;

class Tst34RunAction : public G4UserRunAction
{
  public:
  
  
  Tst34RunAction();
  ~Tst34RunAction();
  void BeginOfRunAction(const G4Run* aRun);
  void EndOfRunAction(const G4Run* aRun);
  void setRunID (int i) {runID=i;};
  private:
    	G4int runID;
};


#endif
