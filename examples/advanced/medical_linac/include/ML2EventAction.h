#ifndef CML2EventActionH
#define CML2EventActionH

#include "G4UserEventAction.hh"
#include "globals.hh"

class CML2EventAction: public G4UserEventAction
{
public:
	CML2EventAction(void);
	~CML2EventAction(void);
  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

private:
  G4String drawFlag; //Visualisation flag
};


#endif
