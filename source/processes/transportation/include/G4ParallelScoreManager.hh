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
// $Id: G4ParallelScoreManager.hh,v 1.3 2002-04-09 17:40:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelScoreManager
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelScoreManager_hh
#define G4ParallelScoreManager_hh

#include "globals.hh"

class G4VPhysicalVolume;
class G4ParallelManager;
class G4VPScorer;
class G4PScoreProcess;

class G4ParallelScoreManager
{

public:  // with description

  G4ParallelScoreManager(G4VPhysicalVolume &worldvolume,
			 const G4String &particlename,
			 G4VPScorer &scorer);
  ~G4ParallelScoreManager();

  G4PScoreProcess *CreateParallelScoreProcess();
  void Initialize();

private:

  G4ParallelScoreManager(const G4ParallelScoreManager &);
  G4ParallelScoreManager &operator=(const G4ParallelScoreManager&);

private:

  G4ParallelManager &fParallelManager;
  G4VPScorer &fPScorer;
  G4PScoreProcess *fPScorerProcess;
};

#endif
