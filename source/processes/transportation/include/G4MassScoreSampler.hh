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
// $Id: G4MassScoreSampler.hh,v 1.1 2002-05-31 10:16:01 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MassScoreSampler
//
// Class description:
//
// A user should use this class to set up scoring in the "mass" geometry.
// The user must create an object of this kind and initialise it.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassScoreSampler_hh
#define G4MassScoreSampler_hh G4MassScoreSampler_hh

#include "globals.hh"
#include "G4VSampler.hh"

class G4VProcess;
class G4VPScorer;
class G4MScoreProcess;

class G4MassScoreSampler : public G4VSampler
{

public:  // with description

  G4MassScoreSampler(G4VPScorer &ascorer, const G4String &particlename);
    // a G4MassScoreSampler for a particle type

  ~G4MassScoreSampler();
    // delete G4MScoreProcess if constructed

  G4MScoreProcess *CreateMassScoreProcess();
    // create the mass score process 
    // don't use it if you use Initialize()

  void Initialize();
    // the G4MassScoreSampler has to be initialised after
    // the initialisation of the G4RunManager !

private:

  G4MassScoreSampler(const G4MassScoreSampler &);
  G4MassScoreSampler &operator=(const G4MassScoreSampler &);

private:

  G4VPScorer &fScorer;
  G4String fParticleName;
  G4MScoreProcess *fMScoreProcess;
};

#endif
