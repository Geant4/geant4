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
// $Id: G4MassScoreManager.hh,v 1.2 2002-04-09 17:40:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MassScoreManager
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassScoreManager_hh
#define G4MassScoreManager_hh

#include "globals.hh"

class G4VProcess;
class G4VPScorer;
class G4MScoreProcess;

class G4MassScoreManager
{

public:  // with description

  G4MassScoreManager(G4VPScorer &ascorer, const G4String &particlename);
  ~G4MassScoreManager();

  G4MScoreProcess *CreateMassScoreProcess();
  void Initialize();

private:

  G4MassScoreManager(const G4MassScoreManager &);
  G4MassScoreManager &operator=(const G4MassScoreManager &);

private:

  G4VPScorer &fScorer;
  G4String fParticleName;
  G4MScoreProcess *fMScoreProcess;
};

#endif
