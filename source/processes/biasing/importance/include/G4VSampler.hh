//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4VSampler.hh 66241 2012-12-13 18:34:42Z gunter $
//
// ----------------------------------------------------------------------
// Class G4VSampler
//
// Class description:
//
// This interface describes a configurable sampler.
// It applies to a given particle type.
// Concrete classes with this interface may be used for 
// scoring, importance sampling and weight cutoff (weight roulette).

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VSampler_hh
#define G4VSampler_hh G4VSampler_hh

#include "G4Types.hh"
#include "G4PlaceOfAction.hh"

class G4VPhysicalVolume;
class G4VImportanceAlgorithm;
class G4VIStore;
class G4VWeightWindowAlgorithm;
class G4VWeightWindowStore;
//class G4VScorer;

class G4VSampler
{

public:  // with description
  
  G4VSampler();
  virtual ~G4VSampler();

  //  virtual void PrepareScoring(G4VScorer *Scorer) = 0;
  //  virtual void PrepareScoring() = 0;

  virtual void PrepareImportanceSampling(G4VIStore *istore,
                                         const G4VImportanceAlgorithm 
                                         *ialg = 0) = 0;


  virtual void PrepareWeightRoulett(G4double wsurvive = 0.5, 
                                    G4double wlimit = 0.25,
                                    G4double isource = 1) = 0;

  virtual void PrepareWeightWindow(G4VWeightWindowStore *wwstore,
                                   G4VWeightWindowAlgorithm *wwAlg = 0,
                                   G4PlaceOfAction placeOfAction = 
                                   onBoundary) = 0;

  virtual void Configure() = 0;

  virtual void ClearSampling() = 0;
    // clear the sampler and remove the processes

  virtual G4bool IsConfigured() const = 0;
    // check if some initialization hase already been done

  virtual void SetParallel(G4bool paraflag) = 0;

};
  
#endif
