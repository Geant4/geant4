//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ParallelImportanceSampler.hh,v 1.3 2002-09-02 13:27:26 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelImportanceSampler
//
// Class description:
//
// A user should use this class to set up importance sampling
// in a "parallel" geometry.
// The user must create an object of this kind and initialise it.


// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelImportanceSampler_hh
#define G4ParallelImportanceSampler_hh G4ParallelImportanceSampler_hh 

#include "globals.hh"
#include "G4VSampler.hh"

class G4ParallelWorld;
class G4VIStore;
class G4VImportanceAlgorithm;
class G4VImportanceSplitExaminer;
class G4ParallelImportanceProcess;
class G4VTrackTerminator;
class G4VPhysicalVolume;

class G4ParallelImportanceSampler : public G4VSampler
{

public:  // with description

  G4ParallelImportanceSampler(G4VPhysicalVolume &worldvolume,
			      G4VIStore &is, 
			      const G4String &particlename,
			      const G4VImportanceAlgorithm *ialg = 0);
    // if *ialg = 0: use the G4ImportanceAlgorithm 
    // use a customised  importance algorithm derived from
    // G4VImportanceAlgorithm  
    // create G4ParallelWorld and G4ImportanceSplitExaminer

public: // used internally by importance sampling 

  G4ParallelImportanceSampler(G4VIStore &is, 
			      const G4String &particlename,
			      G4ParallelWorld &pworld,
			      const G4VImportanceAlgorithm *ialg = 0);
  
  ~G4ParallelImportanceSampler();
    // delete G4ParallelWorld and G4ImportanceSplitExaminer
    // and G4ImportanceAgorithm and G4ParallelImportanceProcess
    // if created

public:  // with description

  G4ParallelImportanceProcess *CreateParallelImportanceProcess();
    // create the parallel importance process 
    // don't use it if you use Initialize()

  void Initialize();
    // the G4ParallelImportanceSampler has to be initialised after
    // the initialisation of the G4RunManager !
  
  void SetTrackTerminator(G4VTrackTerminator *tt){
    fTrackTerminator = tt;
  }
  
 private:

  G4ParallelImportanceSampler(const G4ParallelImportanceSampler &);
  G4ParallelImportanceSampler &operator=(const G4ParallelImportanceSampler &);

private:

  G4String fParticleName;
  G4ParallelWorld *fParallelWorld;
  G4VTrackTerminator *fTrackTerminator;
  G4bool fCreatedPW;

  G4bool fDeleteAlg;
  const G4VImportanceAlgorithm *fIalgorithm;

  G4VImportanceSplitExaminer *fExaminer;
  G4ParallelImportanceProcess *fParallelImportanceProcess;
};

#endif
