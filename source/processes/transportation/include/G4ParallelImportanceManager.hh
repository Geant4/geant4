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
// $Id: G4ParallelImportanceManager.hh,v 1.4 2002-04-09 17:40:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelImportanceManager
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelImportanceManager_hh
#define G4ParallelImportanceManager_hh

#include "globals.hh"

class G4ParallelManager;
class G4VIStore;
class G4VImportanceAlgorithm;
class G4VImportanceSampler;
class G4ParallelImportanceProcess;

class G4ParallelImportanceManager
{

public:  // with description

  G4ParallelImportanceManager(G4VIStore &is, 
			      const G4String &particlename);
  G4ParallelImportanceManager(G4VIStore &is, 
			      G4ParallelManager &pmanager);
  G4ParallelImportanceManager(G4VIStore &is, 
			      const G4String &particlename,
			      G4VImportanceAlgorithm &ialg);
  G4ParallelImportanceManager(G4VIStore &is, 
			      G4VImportanceAlgorithm &ialg,
			      G4ParallelManager &pmanager);
  
  virtual ~G4ParallelImportanceManager();

  G4ParallelImportanceProcess *CreateParallelImportanceProcess();
  void Initialize();

private:

  G4ParallelImportanceManager(const G4ParallelImportanceManager &);
  G4ParallelImportanceManager &operator=(const G4ParallelImportanceManager &);

private:

  G4ParallelManager &fParallelManager;
  bool fCreatedPM;
  G4VImportanceAlgorithm &fIalgorithm;
  G4bool fDeleteAlg;

  G4VImportanceSampler *fSampler;
  G4ParallelImportanceProcess *fParallelImportanceProcess;
};

#endif
