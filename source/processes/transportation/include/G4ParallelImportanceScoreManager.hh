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
// $Id: G4ParallelImportanceScoreManager.hh,v 1.3 2002-04-09 17:40:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelImportanceScoreManager
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelImportanceScoreManager_hh
#define G4ParallelImportanceScoreManager_hh

#include "globals.hh"

class G4VIStore;
class G4VPScorer;
class G4VImportanceAlgorithm;
class G4PScoreProcess;
class G4ParallelManager;
class G4ParallelImportanceManager;

class G4ParallelImportanceScoreManager
{

public:  // with description

  G4ParallelImportanceScoreManager(G4VIStore &iw, 
				   G4VPScorer &ascorer,
				   const G4String &particlename);
  G4ParallelImportanceScoreManager(G4VIStore &iw, 
				   G4VPScorer &ascorer,
				   const G4String &particlename,
				   G4VImportanceAlgorithm &ialg);
  ~G4ParallelImportanceScoreManager();

  G4PScoreProcess *CreateParallelScoreProcess();
  void Initialize();

private:

  G4ParallelImportanceScoreManager(const G4ParallelImportanceScoreManager &);
  G4ParallelImportanceScoreManager &
  operator=(const G4ParallelImportanceScoreManager &);

private:

  G4ParallelManager &fParallelManager;
  G4ParallelImportanceManager &fParallelImportanceManager;
  G4VPScorer &fPScorer;
  G4PScoreProcess *fPScoreProcess;
};

#endif
