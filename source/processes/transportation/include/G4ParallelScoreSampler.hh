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
// $Id: G4ParallelScoreSampler.hh,v 1.3 2002-09-02 13:27:26 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelScoreSampler
//
// Class description:
//
// A user should use this class to set up scoring in a "parallel" 
// geometry.
// The user must create an object of this kind and initialise it.


// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelScoreSampler_hh
#define G4ParallelScoreSampler_hh G4ParallelScoreSampler_hh

#include "globals.hh"
#include "G4VSampler.hh"

class G4ParallelWorld;
class G4VPhysicalVolume;
class G4ParallelTransport;
class G4VPScorer;
class G4PScoreProcess;

class G4ParallelScoreSampler : public G4VSampler
{

public:  // with description

  G4ParallelScoreSampler(G4VPhysicalVolume &worldvolume,
			 const G4String &particlename,
			 G4VPScorer &scorer);
    // create and initalise members 

  ~G4ParallelScoreSampler();
    // delete constructed members

  G4PScoreProcess *CreateParallelScoreProcess();
    // create the parallel score process 
    // don't use it if you use Initialize()
  G4ParallelTransport *CreateParallelTransport();
    // get the G4ParallelTransport
 
  void Initialize();
    // the G4ParallelScoreSampler has to be initialised after
    // the initialisation of the G4RunManager !


private:

  G4ParallelScoreSampler(const G4ParallelScoreSampler &);
  G4ParallelScoreSampler &operator=(const G4ParallelScoreSampler&);

private:
  G4String fParticleName;
  G4ParallelWorld *fParallelWorld;
  G4VPScorer &fPScorer;
  G4PScoreProcess *fPScorerProcess;
  G4ParallelTransport *fParallelTransport;
};

#endif
