// ----------------------------------------------------------------------
// Class G4PImportanceWWindowScoreManager
//
// Class description:
//


// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4PImportanceWWindowScoreManager_hh
#define G4PImportanceWWindowScoreManager_hh G4PImportanceWWindowScoreManager_hh

#include "globals.hh"

class G4VIStore;
class G4VPScorer;
class G4VImportanceAlgorithm;
class G4PScoreProcess;
class G4ParallelManager;
class G4ParallelImportanceManager;
class G4VProcess;
class G4VWeightWindowAlgorithm;

class G4PImportanceWWindowScoreManager
{

public:  // with description

  G4PImportanceWWindowScoreManager(G4VIStore &iw, 
				   G4VPScorer &ascorer,
				   const G4String &particlename,
				   G4VWeightWindowAlgorithm &wwalg);
    // use G4ImportanceAlgorithm, construct and initalise
    // used objects
 
  G4PImportanceWWindowScoreManager(G4VIStore &iw, 
				   G4VPScorer &ascorer,
				   const G4String &particlename,
				   G4VWeightWindowAlgorithm &wwalg,
				   G4VImportanceAlgorithm &ialg);
    // use a customised  importance algorithm derived from
    // G4VImportanceAlgorithm
  

  ~G4PImportanceWWindowScoreManager();
    // delete constructed objects

  G4PScoreProcess *CreateParallelScoreProcess();
    // create the parallel score process 
    // don't use it if you use Initialize()

  G4VProcess *CreateWeightWindowProcess();

  void Initialize();
    // the G4MassImportanceScoreManager has to be initialised after
    // the initialisation of the G4RunManager !

private:

  G4PImportanceWWindowScoreManager(const G4PImportanceWWindowScoreManager &);
  G4PImportanceWWindowScoreManager &
  operator=(const G4PImportanceWWindowScoreManager &);

private:

  G4ParallelManager &fParallelManager;
  G4ParallelImportanceManager &fParallelImportanceManager;
  G4VPScorer &fPScorer;
  G4VIStore &fIstore;
  G4PScoreProcess *fPScoreProcess;
  G4VProcess *fPWeightWindowProcess;
  G4VWeightWindowAlgorithm &fWWAlgorithm;
};

#endif
