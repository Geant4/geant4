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
// $Id: G4PImportanceWWindowScoreManager.hh,v 1.4 2002-05-30 15:57:26 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PImportanceWWindowScoreManager
//
// Class description:
//
// This class manages the importance sampling and importance weight window
// sampling together with scoring in a parallel geometry.
// The user must create an object of this kind and initialise it.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4PImportanceWWindowScoreManager_hh
#define G4PImportanceWWindowScoreManager_hh G4PImportanceWWindowScoreManager_hh

#include "globals.hh"
#include "G4VSampler.hh"

class G4VIStore;
class G4VPScorer;
class G4VImportanceAlgorithm;
class G4PScoreProcess;
class G4ParallelManager;
class G4ParallelImportanceManager;
class G4VProcess;
class G4VWeightWindowAlgorithm;

class G4PImportanceWWindowScoreManager : 
  public G4VSampler
{

public:  // with description

  G4PImportanceWWindowScoreManager(G4VIStore &iw, 
				   G4VPScorer &ascorer,
				   const G4String &particlename,
				   G4VWeightWindowAlgorithm &wwalg,
				   const G4VImportanceAlgorithm *ialg = 0);
    // if *ialg = 0: use G4ImportanceAlgorithm
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
