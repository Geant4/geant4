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
// $Id: G4ParallelImportanceScoreSampler.hh,v 1.1 2002-05-31 10:16:01 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelImportanceScoreSampler
//
// Class description:
//
// A user should use this class to set up importance sampling and scoring
// in a "parallel" geometry.
// The user must create an object of this kind and initialise it.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelImportanceScoreSampler_hh
#define G4ParallelImportanceScoreSampler_hh G4ParallelImportanceScoreSampler_hh

#include "globals.hh"
#include "G4VSampler.hh"

class G4VIStore;
class G4VPScorer;
class G4VImportanceAlgorithm;
class G4PScoreProcess;
class G4ParallelWorld;
class G4ParallelImportanceSampler;

class G4ParallelImportanceScoreSampler : 
  public G4VSampler
{

public:  // with description

 
  G4ParallelImportanceScoreSampler(G4VIStore &iw, 
				   G4VPScorer &ascorer,
				   const G4String &particlename,
				   const G4VImportanceAlgorithm *ialg = 0);
    // if *ialg = 0: use G4ImportanceAlgorithm,
    // use a customised  importance algorithm derived from
    // G4VImportanceAlgorithm
  

  ~G4ParallelImportanceScoreSampler();
    // delete constructed objects

  G4PScoreProcess *CreateParallelScoreProcess();
    // create the parallel score process 
    // don't use it if you use Initialize()

  void Initialize();
    // the G4MassImportanceScoreSampler has to be initialised after
    // the initialisation of the G4RunManager !

private:

  G4ParallelImportanceScoreSampler(const G4ParallelImportanceScoreSampler &);
  G4ParallelImportanceScoreSampler &
  operator=(const G4ParallelImportanceScoreSampler &);

private:
  G4String fParticleName;
  G4ParallelWorld &fParallelWorld;
  G4ParallelImportanceSampler &fParallelImportanceSampler;
  G4VPScorer &fPScorer;
  G4PScoreProcess *fPScoreProcess;
};

#endif
