// ----------------------------------------------------------------------
// Class G4ParallelWeightWindowProcess
//
// Class description:
//


// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelWeightWindowProcess_hh
#define G4ParallelWeightWindowProcess_hh G4ParallelWeightWindowProcess_hh

#include "g4std/strstream"

#include "G4VProcess.hh"
#include "globals.hh"
#include "G4Nsplit_Weight.hh"
#include "G4ImportancePostStepDoIt.hh"

class G4VIStore;;
class G4VParallelStepper;
class G4VWeightWindowAlgorithm;

class G4ParallelWeightWindowProcess : public G4VProcess
{

public:  // with description

  G4ParallelWeightWindowProcess(G4VIStore &aIstore,
				G4VParallelStepper &aStepper,
				G4VWeightWindowAlgorithm &aWWAlgorithm,
				const G4String &aName = "ParallelWeightWindowProcess");
    // create G4ParticleChange

  virtual ~G4ParallelWeightWindowProcess();
    // delete G4ParticleChange



  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double   previousStepSize,
				       G4ForceCondition* condition);


  virtual G4VParticleChange *PostStepDoIt(const G4Track&, const G4Step&);


  void StartTracking();
  void EndTracking();

public:  // without description
  
  //  no operation in  AtRestDoIt and  AlongStepDoIt

  G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
                                        G4double  ,
                                        G4double  ,
                                        G4double& ,
                                        G4GPILSelection*) {return -1.0;}
  G4double AtRestGetPhysicalInteractionLength(const G4Track& ,
                                        G4ForceCondition*) {return -1.0;}
  
  G4VParticleChange*  AtRestDoIt(const G4Track&, const G4Step&) {return 0;}
  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) {return 0;}

protected:

  virtual void Error(const G4String &m);
  virtual void Warning(const G4String &m);

  G4ParticleChange *fParticleChange;


private:

  G4ParallelWeightWindowProcess(const G4ParallelWeightWindowProcess &);
  G4ParallelWeightWindowProcess &operator=(const G4ParallelWeightWindowProcess &);

private:

  G4VIStore &fIStore;
  G4VParallelStepper &fPStepper;
  G4VWeightWindowAlgorithm &fWWAlgorithm;
  G4Nsplit_Weight fNsplit_Weight;
  G4ImportancePostStepDoIt fImportancePostStepDoIt;
  G4bool fInitStep;
};

#endif
