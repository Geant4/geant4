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
// $Id: G4ParallelImportanceProcess.hh,v 1.3 2002-04-10 13:14:17 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelImportanceProcess
//
// Class description:
//
// Usd internally by importance sampling in a "parallel" geometry.
// This is a G4ParallelTransport that also does importance
// sampling in the "parallel" geometry.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelImportanceProcess_hh
#define G4ParallelImportanceProcess_hh G4ParallelImportanceProcess_hh

#include "G4ParallelTransport.hh"
#include "G4ImportancePostStepDoIt.hh"

class G4VImportanceSampler;
class G4Nsplit_Weight;

class G4ParallelImportanceProcess : public G4ParallelTransport
{

public:  // with description

  G4ParallelImportanceProcess(const G4VImportanceSampler &aImportanceSampler,
			      G4VPGeoDriver &pgeodriver, 
			      G4VParallelStepper &aStepper,
			      const G4String &aName = "ParallelImportanceProcess");  
    // initialise G4ParallelTransport and members

  G4VParticleChange *PostStepDoIt(const G4Track&,
				  const G4Step&);
    // do the "parallel transport" and importance sampling.

private:

  G4ParallelImportanceProcess(const G4ParallelImportanceProcess &);
  G4ParallelImportanceProcess &operator=(const G4ParallelImportanceProcess &);
  
  void Error(const G4String &m);

private:

  G4ParticleChange *fParticleChange;
  const G4VImportanceSampler &fImportanceSampler;  
  G4ImportancePostStepDoIt fImportancePostStepDoIt;
};

#endif
