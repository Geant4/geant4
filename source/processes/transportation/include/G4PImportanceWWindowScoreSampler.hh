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
// $Id: G4PImportanceWWindowScoreSampler.hh,v 1.3 2002-09-02 13:27:26 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PImportanceWWindowScoreSampler
//
// Class description:
//
// This class manages the importance sampling and importance weight window
// sampling together with scoring in a parallel geometry.
// The user must create an object of this kind and initialise it.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4PImportanceWWindowScoreSampler_hh
#define G4PImportanceWWindowScoreSampler_hh G4PImportanceWWindowScoreSampler_hh

#include "globals.hh"
#include "G4VSampler.hh"

class G4VIStore;
class G4VPScorer;
class G4VImportanceAlgorithm;
class G4PScoreProcess;
class G4ParallelWorld;
class G4ParallelImportanceSampler;
class G4VProcess;
class G4VWeightWindowAlgorithm;

class G4PImportanceWWindowScoreSampler : 
  public G4VSampler
{

public:  // with description

  G4PImportanceWWindowScoreSampler(G4VPhysicalVolume &worldvolume,
				   G4VIStore &iw, 
				   G4VPScorer &ascorer,
				   const G4String &particlename,
				   G4VWeightWindowAlgorithm &wwalg,
				   const G4VImportanceAlgorithm *ialg = 0);
    // if *ialg = 0: use G4ImportanceAlgorithm
    // use a customised  importance algorithm derived from
    // G4VImportanceAlgorithm
  

  ~G4PImportanceWWindowScoreSampler();
    // delete constructed objects

  G4PScoreProcess *CreateParallelScoreProcess();
    // create the parallel score process 
    // don't use it if you use Initialize()

  G4VProcess *CreateWeightWindowProcess();

  void Initialize();
    // the G4MassImportanceScoreSampler has to be initialised after
    // the initialisation of the G4RunManager !

private:

  G4PImportanceWWindowScoreSampler(const G4PImportanceWWindowScoreSampler &);
  G4PImportanceWWindowScoreSampler &
  operator=(const G4PImportanceWWindowScoreSampler &);

private:
  G4String fParticleName;
  G4ParallelWorld *fParallelWorld;
  G4ParallelImportanceSampler *fParallelImportanceSampler;
  G4VPScorer &fPScorer;
  G4VIStore &fIstore;
  G4PScoreProcess *fPScoreProcess;
  G4VProcess *fPWeightWindowProcess;
  G4VWeightWindowAlgorithm &fWWAlgorithm;
};

#endif
