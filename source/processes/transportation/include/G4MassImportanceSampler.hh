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
// $Id: G4MassImportanceSampler.hh,v 1.1 2002-05-31 10:16:01 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MassImportanceSampler
//
// Class description:
//
// A user should use this class to set up importance sampling
// in the "mass" geometry.
// The user must create an object of this kind and initialise it.


// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassImportanceSampler_hh
#define G4MassImportanceSampler_hh G4MassImportanceSampler_hh

#include "globals.hh"
#include "G4VSampler.hh"

class G4VIStore;
class G4MassImportanceProcess;
class G4VImportanceAlgorithm;
class G4VProcess;

class G4MassImportanceSampler : public G4VSampler
{

public:  // with description

  G4MassImportanceSampler(G4VIStore &aIstore,
			  const G4String &particlename,
			  const G4VImportanceAlgorithm *algorithm = 0);
    // if *algorithm = 0: use the G4ImportanceAlgorithm 
    // else: use a customised  importance algorithm derived from
    // G4VImportanceAlgorithm

  ~G4MassImportanceSampler();
    // delete G4MassImportanceProcess and G4VImportanceAlgorithm
    // if constructed by this object

  G4MassImportanceProcess *CreateMassImportanceProcess();
    // create the mass importance process 
    // don't use it if you use Initialize()
  
  void Initialize();
    // the G4MassImportanceSampler has to be initialised after
    // the initialisation of the G4RunManager !
  

private:

  G4MassImportanceSampler(const G4MassImportanceSampler &);
  G4MassImportanceSampler &operator=(const G4MassImportanceSampler &);

private:

  G4VIStore &fIStore;
  G4String fParticleName;
  G4MassImportanceProcess *fMassImportanceProcess;
  G4bool fCreatedAlgorithm;
  const G4VImportanceAlgorithm *fAlgorithm;
};

#endif
