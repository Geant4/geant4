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
// $Id: G4MassScoreManager.hh,v 1.4 2002-05-30 11:14:38 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MassScoreManager
//
// Class description:
//
// A user should use this class to set up scoring in the "mass" geometry.
// Create an object and initialise it.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassScoreManager_hh
#define G4MassScoreManager_hh G4MassScoreManager_hh

#include "globals.hh"
#include "G4VImportanceScoreConstructor.hh"

class G4VProcess;
class G4VPScorer;
class G4MScoreProcess;

class G4MassScoreManager : public G4VImportanceScoreConstructor
{

public:  // with description

  G4MassScoreManager(G4VPScorer &ascorer, const G4String &particlename);
    // a G4MassScoreManager for a particle type

  ~G4MassScoreManager();
    // delete G4MScoreProcess if constructed

  G4MScoreProcess *CreateMassScoreProcess();
    // create the mass score process 
    // don't use it if you use Initialize()

  void Initialize();
    // the G4MassScoreManager has to be initialised after
    // the initialisation of the G4RunManager !

private:

  G4MassScoreManager(const G4MassScoreManager &);
  G4MassScoreManager &operator=(const G4MassScoreManager &);

private:

  G4VPScorer &fScorer;
  G4String fParticleName;
  G4MScoreProcess *fMScoreProcess;
};

#endif
