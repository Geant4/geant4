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
// $Id: G4ImportanceSampler.hh,v 1.3 2002-04-10 13:13:06 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ImportanceSampler
//
// Class description:
//
// This class is used internally by importance sampling.
// Implementation of a importance sampler (see G4VImportanceSampler).
// This implementation is used for sampling in a "parallel"
// geometry.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ImportanceSampler_hh
#define G4ImportanceSampler_hh G4ImportanceSampler_hh

#include "G4VImportanceSampler.hh"

class G4VImportanceAlgorithm;
class G4VParallelStepper;
class G4ImportanceFinder;
class G4VIStore;

class G4ImportanceSampler: public G4VImportanceSampler
{

public:  // with description

  G4ImportanceSampler(const G4VImportanceAlgorithm &aIalg,
		      const G4VParallelStepper &astepper,
		      const G4VIStore &istore);
    // initialisation and construct G4ImportanceFinder

  ~G4ImportanceSampler();
    // delete G4ImportanceFinder

  G4Nsplit_Weight Sample(G4double w) const; 
    // Get  G4Nsplit_Weight for a given mother track weight.
  

private:

  G4ImportanceSampler(const G4ImportanceSampler &);
  G4ImportanceSampler &operator=(const G4ImportanceSampler &);

private:

  const G4VImportanceAlgorithm &fIalgorithm;
  const G4VParallelStepper &fPStepper;
  const G4ImportanceFinder &fIfinder;
};


#endif
