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
// $Id: G4LowEnergyGammaConversion.hh,v 1.13 2001-09-10 18:05:16 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: A. Forti
//         Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 02 Mar 1999   A. Forti   1st implementation
// 14 Aug 2001   MGP        Major revision according to a design iteration
//
// -------------------------------------------------------------------

// Class description:
// Low Energy process, Photon conversion
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4LOWENERGYGAMMACONVERSION_HH
#define G4LOWENERGYGAMMACONVERSION_HH 1


#include "G4VDiscreteProcess.hh"

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4CrossSectionHandler;
class G4VDataSetAlgorithm;

class G4LowEnergyGammaConversion : public G4VDiscreteProcess {

public:
 
  G4LowEnergyGammaConversion(const G4String& processName ="LowEnConversion");
 
  ~G4LowEnergyGammaConversion();

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
  G4LowEnergyGammaConversion& operator=(const G4LowEnergyGammaConversion &right);
  G4LowEnergyGammaConversion(const G4LowEnergyGammaConversion& );

  G4double ScreenFunction1(G4double screenVariable);
  G4double ScreenFunction2(G4double screenVariable);
  
  G4double lowEnergyLimit;  // low energy limit  applied to the process
  G4double highEnergyLimit; // high energy limit applied to the process

  G4VEMDataSet* meanFreePathTable;
  G4CrossSectionHandler* crossSectionHandler;

  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;

  const G4double smallEnergy;

};

#endif
 








