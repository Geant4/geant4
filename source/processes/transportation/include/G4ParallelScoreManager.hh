#ifndef G4ParallelScoreManager_hh
#define G4ParallelScoreManager_hh G4ParallelScoreManager_hh

#include "G4ParallelManager.hh"


class G4VPScorer;
class G4PScoreProcess;

class G4ParallelScoreManager : private G4ParallelManager{
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
  G4VPScorer &fPScorer;
  G4PScoreProcess *fPScorerProcess;
};

#endif
