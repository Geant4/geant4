#ifndef G4ParallelScoreManager_hh
#define G4ParallelScoreManager_hh G4ParallelScoreManager_hh

#include "globals.hh"

class G4VPhysicalVolume;
class G4ParallelManager;
class G4VPScorer;
class G4PScoreProcess;

class G4ParallelScoreManager{
public:
  G4ParallelScoreManager(G4VPhysicalVolume &worldvolume,
			 const G4String &particlename,
			 G4VPScorer &scorer);

  ~G4ParallelScoreManager();
  G4PScoreProcess *CreateParallelScoreProcess();
  void Initialize();

private:
  G4ParallelScoreManager(const G4ParallelScoreManager &);
  G4ParallelScoreManager &operator=(const G4ParallelScoreManager&);
  G4ParallelManager &fParallelManager;
  G4VPScorer &fPScorer;
  G4PScoreProcess *fPScorerProcess;
};

#endif
