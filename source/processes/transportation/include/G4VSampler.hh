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
// $Id: G4VSampler.hh,v 1.6 2002-10-22 13:25:57 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VSampler
//
// Class description:
//
// This interface discribes a configurable sampler.
// It applies to a given particle type.
// Concrete classes with this interface may be used for 
// scoring, importance sampling and weigth cutoff (weight roulett).
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VSampler_hh
#define G4VSampler_hh G4VSampler_hh


#include "globals.hh"

class G4VPhysicalVolume;
class G4VImportanceAlgorithm;
class G4VIStore;
class G4VScorer;


class G4VSampler {

public:  
  
  G4VSampler();
  virtual ~G4VSampler();

  virtual void PrepareScoring(G4VScorer *Scorer) = 0;

  virtual void PrepareImportanceSampling(G4VIStore *istore,
					 const G4VImportanceAlgorithm 
					 *ialg = 0) = 0;


  virtual void PrepareWeightRoulett(G4double wsurvive = 0.5, 
				    G4double wlimit = 0.25,
				    G4double isource = 1) = 0;

  virtual void Configure() = 0;

  virtual void ClearSampling() = 0;
    // clear the sampler and remove the processes

  virtual G4bool IsConfigured() const = 0;
    // check if some initialization hase already been done
};
  
#endif

