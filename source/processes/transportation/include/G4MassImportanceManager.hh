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
// $Id: G4MassImportanceManager.hh,v 1.2 2002-04-09 17:40:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MassImportanceManager
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassImportanceManager_hh
#define G4MassImportanceManager_hh

#include "globals.hh"

class G4VIStore;
class G4MassImportanceProcess;
class G4VImportanceAlgorithm;
class G4VProcess;

class G4MassImportanceManager
{

public:  // with description

  G4MassImportanceManager(G4VIStore &aIstore,
			  const G4String &particlename);

  G4MassImportanceManager(G4VIStore &aIstore,
			  const G4String &particlename,
			  const G4VImportanceAlgorithm &algorithm);
  
  ~G4MassImportanceManager();

  G4MassImportanceProcess *CreateMassImportanceProcess();
  void Initialize();

private:

  G4MassImportanceManager(const G4MassImportanceManager &);
  G4MassImportanceManager &operator=(const G4MassImportanceManager &);

private:

  G4VIStore &fIStore;
  G4String fParticleName;
  const G4VImportanceAlgorithm &fAlgorithm;
  bool fCreatedAlgorithm;
  G4MassImportanceProcess *fMassImportanceProcess;
};

#endif
