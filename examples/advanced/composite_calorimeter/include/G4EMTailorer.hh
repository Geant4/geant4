#ifndef G4EMTailorer_h
#define G4EMTailorer_h

class G4EMBuilder;
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

class G4EMTailorer: public G4UImessenger
{
  public:
    G4EMTailorer(G4EMBuilder * af);
    
    virtual ~G4EMTailorer()
    {
      delete theSynch;
      delete theGN;
    }

    void SetNewValue(G4UIcommand* aComm, G4String aS);
    
  private:
    G4EMBuilder * theB;
    G4UIcmdWithAString * theSynch;
    G4UIcmdWithAString * theGN;
    G4UIdirectory *aDir1;
    G4UIdirectory *aDir2;
};

#endif
