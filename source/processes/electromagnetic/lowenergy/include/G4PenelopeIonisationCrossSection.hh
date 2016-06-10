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
// $Id: G4PenelopeIonisationCrossSection.hh 66241 2012-12-13 18:34:42Z gunter $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 14 Mar 2012   L. Pandola   1st implementation. 
//
// -------------------------------------------------------------------
//
//! Class description:
//!  Access to Penelope ionisation cross sections for e- and e+. To be 
//!  registered as a specific model in G4UAtomicDeexcitation.
//!
//! NOTICE: working only for e- at the moment (no interface available for
//!         e+)
//!
// -------------------------------------------------------------------

#ifndef G4PENELOPEIONISATIONCROSSSECTION_HH
#define G4PENELOPEIONISATIONCROSSSECTION_HH 1

#include "globals.hh"
#include "G4VhShellCrossSection.hh"
#include <map>

class G4AtomicTransitionManager;
class G4PenelopeOscillatorManager;
class G4PenelopeOscillator;
class G4PenelopeCrossSection;
class G4PenelopeIonisationXSHandler;

class G4PenelopeIonisationCrossSection : public G4VhShellCrossSection
{
public:  
  //! Constructor. 
  G4PenelopeIonisationCrossSection();
  
  //! Destructor. Clean all tables.
  ~G4PenelopeIonisationCrossSection();

  //! Purely virtual method from the base interface. Returns the cross 
  //! section for all levels of element Z in material mat at the 
  //! given energy
  std::vector<G4double> GetCrossSection(G4int Z,
					G4double incidentEnergy,
					G4double mass,
					G4double deltaEnergy,
					const G4Material* mat);
  
  //! Purely virtual method from the base interface. Returns the 
  //! cross section for the given shell in the element Z of material 
  //! mat at the specified energy
  G4double CrossSection(G4int Z,
			G4AtomicShellEnumerator shell,
			G4double incidentEnergy,
			G4double mass,
			const G4Material* mat);

  //! Purely virtual method from the base interface. Returns the 
  //! shell ionisation probabilities for the given Z in the
  //! material mat at the specified energy.
  std::vector<G4double> Probabilities(G4int Z,
				      G4double incidentEnergy,
				      G4double mass,
				      G4double deltaEnergy,
				      const G4Material* mat) ;
  //! Getter/setter for the verbosity level
  void SetVerbosityLevel(G4int vl){verboseLevel = vl;};
  G4int GetVerbosityLevel(){return verboseLevel;};

private:  
  G4PenelopeIonisationCrossSection & operator=(const G4PenelopeIonisationCrossSection &right);
  G4PenelopeIonisationCrossSection(const G4PenelopeIonisationCrossSection&);

  //Oscillator manager
  G4PenelopeOscillatorManager* oscManager;

  G4int verboseLevel;

  //!The shells in Penelope are organized per *material*, rather than per 
  //!element, so given a material one has to find the proper index for the 
  //!given Z and shellID. An appropriate look-up table is used to avoid 
  //!recalculation.
  G4int FindShellIDIndex(const G4Material* mat,G4int Z,G4AtomicShellEnumerator shell);
  std::map< std::pair<const G4Material*,G4int>, G4DataVector*> *shellIDTable;

  G4int nMaxLevels;

  G4double fLowEnergyLimit;
  G4double fHighEnergyLimit;

  G4PenelopeIonisationXSHandler* theCrossSectionHandler;
  const G4AtomicTransitionManager* transitionManager;
};

#endif

