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
// $Id: G4PenelopeGammaConversion.hh,v 1.1 2002-12-06 16:25:00 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: L. Pandola
//         
// History:
// -----------
// 02 Dec 2002   L.Pandola  1st implementation
// -------------------------------------------------------------------

// Class description:
// Low Energy process, Photon conversion
// Penelope model

// -------------------------------------------------------------------

#ifndef G4PENELOPEGAMMACONVERSION_HH
#define G4PENELOPEGAMMACONVERSION_HH 1

#include "G4VDiscreteProcess.hh"

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4VCrossSectionHandler;
class G4VRangeTest;

class G4PenelopeGammaConversion : public G4VDiscreteProcess {

public:
 
  G4PenelopeGammaConversion(const G4String& processName ="PenConversion");
 
  ~G4PenelopeGammaConversion();

  G4bool IsApplicable(const G4ParticleDefinition& photon);

  void BuildPhysicsTable(const G4ParticleDefinition& photon);


  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
 
  // For testing purpose only
  G4double DumpMeanFreePath(const G4Track& aTrack, 
			    G4double previousStepSize, 
			    G4ForceCondition* condition) 
  { return GetMeanFreePath(aTrack, previousStepSize, condition); }

protected:
  
  G4double GetMeanFreePath(const G4Track& aTrack, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);

private:

  // Hide copy constructor and assignment operator as private 
  G4PenelopeGammaConversion& operator=(const G4PenelopeGammaConversion &right);
  G4PenelopeGammaConversion(const G4PenelopeGammaConversion& );

  G4double GetScreeningRadius(G4double Z); 
  G4double ScreenFunction(G4double screenVariable,G4int icase); 
  G4double CoulombCorrection(G4double ZAlpha); 
  G4double LowEnergyCorrection(G4double ZAlpha,G4double eki); 
  
  G4double lowEnergyLimit;  // low energy limit  applied to the process
  G4double highEnergyLimit; // high energy limit applied to the process

  G4VEMDataSet* meanFreePathTable;
  G4VCrossSectionHandler* crossSectionHandler;

  G4VRangeTest* rangeTest;

  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;

  const G4double smallEnergy;

};

#endif
 








