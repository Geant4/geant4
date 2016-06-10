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
// $Id: G4LivermoreIonisationCrossSection.hh 66241 2012-12-13 18:34:42Z gunter $
//
//
// Author: Vladimir Ivanchenko
//
// History:
// --------
// 31 May 2011 V.Ivanchenko  The class is created  
// 04 Jul 2011 L Pandola     Comment unused private member
// 09 Mar 2012 L Pandola     Changed signature of methods
// 
// 
// -------------------------------------------------------------------
//
// Class description:
// Access to Livermore ionisation cross sections for e-
// -------------------------------------------------------------------

#ifndef G4LIVERMOREIONISATIONCROSSSECTION_HH
#define G4LIVERMOREIONISATIONCROSSSECTION_HH 1

#include "G4VhShellCrossSection.hh"
#include "globals.hh"
#include "G4AtomicShellEnumerator.hh"
#include <vector>

class G4AtomicTransitionManager;
class G4VCrossSectionHandler;

class G4LivermoreIonisationCrossSection : public G4VhShellCrossSection 
{

public:
  
  G4LivermoreIonisationCrossSection(const G4String& nam = "LivermorePIXE");
  
  virtual ~G4LivermoreIonisationCrossSection();

  void Initialise();

  G4double CrossSection(G4int Z, G4AtomicShellEnumerator shell,
			G4double incidentEnergy,
			G4double mass = 0.0,
			const G4Material* mat = 0);

  std::vector<G4double> GetCrossSection(G4int Z,
					G4double incidentEnergy,
					G4double mass = 0.0,
					G4double deltaEnergy = 0.0,
					const G4Material* mat = 0);

  std::vector<G4double> Probabilities(G4int Z,
				      G4double incidentEnergy,
				      G4double mass = 0.0,
				      G4double deltaEnergy = 0,
				      const G4Material* mat = 0);
    
    
private:
 
  G4LivermoreIonisationCrossSection & operator=(const G4LivermoreIonisationCrossSection &right);
  G4LivermoreIonisationCrossSection(const G4LivermoreIonisationCrossSection&);

  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fLowEnergyLimit;
  G4double fHighEnergyLimit;
 
  //G4bool isInitialised;

  G4int verboseLevel;
 
  G4VCrossSectionHandler* crossSectionHandler;

  const G4AtomicTransitionManager* transitionManager;

};

#endif

