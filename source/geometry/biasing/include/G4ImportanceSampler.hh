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
// $Id: G4ImportanceSampler.hh,v 1.2 2002-04-09 16:23:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ImportanceSampler
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ImportanceSampler_hh
#define G4ImportanceSampler_hh

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
  ~G4ImportanceSampler();

  G4Nsplit_Weight Sample(G4double w) const; 

private:

  G4ImportanceSampler(const G4ImportanceSampler &);
  G4ImportanceSampler &operator=(const G4ImportanceSampler &);

private:

  const G4VImportanceAlgorithm &fIalgorithm;
  const G4VParallelStepper &fPStepper;
  const G4ImportanceFinder &fIfinder;
};


#endif
