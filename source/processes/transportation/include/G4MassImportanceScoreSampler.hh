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
// $Id: G4MassImportanceScoreSampler.hh,v 1.1 2002-05-31 10:16:01 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MassImportanceScoreSampler
//
// Class description:
//
// A user should use this class to set up importance sampling and scoring
// in the "mass" geometry.
// The user must create an object of this kind and initialise it.


// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4MassImportanceScoreSampler_hh
#define G4MassImportanceScoreSampler_hh G4MassImportanceScoreSampler_hh

#include "globals.hh"
#include "G4VSampler.hh"

class G4VIStore;
class G4VPScorer;
class G4VImportanceAlgorithm;
class G4MassImportanceSampler;
class G4MassScoreSampler;

class G4MassImportanceScoreSampler : public G4VSampler
{

public:  // with description

  G4MassImportanceScoreSampler(G4VIStore &aIstore,
			       G4VPScorer &ascorer,
			       const G4String &particlename,
			       const G4VImportanceAlgorithm 
			       *algorithm = 0);
    // if *algorithm = 0: use the G4ImportanceAlgorithm
    // else : use a customised  importance algorithm derived from
    // G4VImportanceAlgorithm  
  

  ~G4MassImportanceScoreSampler();
    // delete G4MassScoreSampler and G4MassImportanceSampler if created

  void Initialize();
    // the G4MassImportanceScoreSampler has to be initialised after
    // the initialisation of the G4RunManager !

private:

  G4MassImportanceScoreSampler(const G4MassImportanceScoreSampler &);
  G4MassImportanceScoreSampler &
  operator=(const G4MassImportanceScoreSampler &);

private:

  G4MassImportanceSampler *fMassImportanceSampler;
  G4MassScoreSampler *fMassScoreSampler;
};

#endif
