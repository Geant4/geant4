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
// $Id: G4QInelastic.hh,v 1.1 2004-03-05 13:26:57 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: M. Kossov (Mikhail.Kossov@cern.ch)
//
// History:
// ----------------------------------------------------------------------------------------
// 20 Dec 2003   M. Kossov  "Hadronic package independent" lepto-nuclear process is created
// ----------------------------------------------------------------------------------------
//
// Class description:
// ==================
// Lepto-nuclear interaction based on the paper:
// M.V.Kossov, Eur.Phys. J. A 14, 377-392 (2002)
// -------------------------------------------------------------------

#ifndef G4LOWENERGYCOMPTON_HH
#define G4LOWENERGYCOMPTON_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4VParticleChange.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4QIsotope.hh"
#include "G4QEnvironment.hh"
#include "G4ParticleTypes.hh"
#include "G4QIsotope.hh"
#include "G4ParticleDefinition.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4QPhotoNuclearCrossSection.hh"

class G4QInelastic : public G4VDiscreteProcess
{
public:
  // Constructors and destructors
  G4QInelastic(const G4String& processName ="CHIPSInelasticReaction");
  virtual ~G4QInelastic();
  // Definition of virtual functions
  G4bool IsApplicable(const G4ParticleDefinition& definition);
  G4double GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition*);
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  // Selector/Modifires
  G4double getLowLimit()      {return lowEnergyLimit;}
  void     setLowLimit(G4double low) {lowEnergyLimit=low;}
private: 
  // Hide copy constructor and assignment operator as private 
  G4QInelastic& operator=(const G4QInelastic& right);
  G4QInelastic(const G4QInelastic& );

  G4double lowEnergyLimit;  // low energy limit applied to the process
  G4QVInelasticCrossSection*   theInelCS;
  G4QPhotoNuclearCrossSection* thePhotonCS;
};

#endif

