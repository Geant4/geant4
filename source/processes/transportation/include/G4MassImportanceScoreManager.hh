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
// $Id: G4MassImportanceScoreManager.hh,v 1.6 2002-05-30 15:57:26 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MassImportanceScoreManager
//
// Class description:
//
// A user should use this class to set up importance sampling and scoring
// in the "mass" geometry.
// The user must create an object of this kind and initialise it.


// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassImportanceScoreManager_hh
#define G4MassImportanceScoreManager_hh G4MassImportanceScoreManager_hh

#include "globals.hh"
#include "G4VSampler.hh"

class G4VIStore;
class G4VPScorer;
class G4VImportanceAlgorithm;
class G4MassImportanceManager;
class G4MassScoreManager;

class G4MassImportanceScoreManager : public G4VSampler
{

public:  // with description

  G4MassImportanceScoreManager(G4VIStore &aIstore,
			       G4VPScorer &ascorer,
			       const G4String &particlename,
			       const G4VImportanceAlgorithm 
			       *algorithm = 0);
    // if *algorithm = 0: use the G4ImportanceAlgorithm
    // else : use a customised  importance algorithm derived from
    // G4VImportanceAlgorithm  
  

  ~G4MassImportanceScoreManager();
    // delete G4MassScoreManager and G4MassImportanceManager if created

  void Initialize();
    // the G4MassImportanceScoreManager has to be initialised after
    // the initialisation of the G4RunManager !

private:

  G4MassImportanceScoreManager(const G4MassImportanceScoreManager &);
  G4MassImportanceScoreManager &
  operator=(const G4MassImportanceScoreManager &);

private:

  G4MassImportanceManager *fMassImportanceManager;
  G4MassScoreManager *fMassScoreManager;
};

#endif
