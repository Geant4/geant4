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
// $Id: G4MassImportanceManager.hh,v 1.3 2002-04-10 13:14:16 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MassImportanceManager
//
// Class description:
//
// A user should use this class to set up importance sampling
// in the "mass" geometry.
// Create an object and initialise it.


// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassImportanceManager_hh
#define G4MassImportanceManager_hh G4MassImportanceManager_hh

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
    // use the G4ImportanceAlgorithm

  G4MassImportanceManager(G4VIStore &aIstore,
			  const G4String &particlename,
			  const G4VImportanceAlgorithm &algorithm);
    // use a customised  importance algorithm derived from
    // G4VImportanceAlgorithm

  ~G4MassImportanceManager();
    // delete G4MassImportanceProcess and G4VImportanceAlgorithm
    // if constructed by this object

  G4MassImportanceProcess *CreateMassImportanceProcess();
    // create the mass importance process 
    // don't use it if you use Initialize()
  
  void Initialize();
    // the G4MassImportanceManager has to be initialised after
    // the initialisation of the G4RunManager !
  

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
