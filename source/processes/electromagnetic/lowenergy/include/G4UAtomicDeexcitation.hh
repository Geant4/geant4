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
// $Id: G4UAtomicDeexcitation.cc,v 1.11 
//
// -------------------------------------------------------------------
//
// Geant4 Header G4UAtomicDeexcitation
//  
// Authors: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// Created 22 April 2010 from old G4AtomicDeexcitation class 
//
// Modified:
// ---------
//  
//
// -------------------------------------------------------------------
//
// Class description:
// Implementation of atomic deexcitation 
//
// -------------------------------------------------------------------

#ifndef G4UAtomicDeexcitation_h
#define G4UAtomicDeexcitation_h 1

#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include <vector>

class G4AtomicTransitionManager;
class G4VhShellCrossSection;
class G4EmCorrections;
class G4Material;

class G4UAtomicDeexcitation : public G4VAtomDeexcitation
{  
public: 
  
  G4UAtomicDeexcitation();
  virtual ~G4UAtomicDeexcitation();
   
  //=================================================================
  // methods that are requested to be implemented by the interface
  //=================================================================

  // initialisation methods
  virtual void InitialiseForNewRun();
  virtual void InitialiseForExtraAtom(G4int Z);


  // Set threshold energy for fluorescence 
  void SetCutForSecondaryPhotons(G4double cut);

  // Set threshold energy for Auger electron production
  void SetCutForAugerElectrons(G4double cut);
  

  // Get atomic shell by shell index, used by discrete processes 
  // (for example, photoelectric), when shell vacancy sampled by the model
  virtual 
  const G4AtomicShell* GetAtomicShell(G4int Z, 
				      G4AtomicShellEnumerator shell);

  // generation of deexcitation for given atom, shell vacancy and cuts
  virtual void GenerateParticles(std::vector<G4DynamicParticle*>* secVect,  
				 const G4AtomicShell*, 
				 G4int Z,
                                 G4double gammaCut,
				 G4double eCut);

  //  access or compute PIXE cross section 
  virtual
  G4double GetShellIonisationCrossSectionPerAtom(const G4ParticleDefinition*, 
						 G4int Z, 
						 G4AtomicShellEnumerator shell,
						 G4double kinE,
                                                 const G4Material* mat = 0);

  //  access or compute PIXE cross section 
  virtual
  G4double ComputeShellIonisationCrossSectionPerAtom(const G4ParticleDefinition*, 
						     G4int Z, 
						     G4AtomicShellEnumerator shell,
						     G4double kinE,
						     const G4Material* mat = 0);

  //=================================================================
  // concrete methods of the deextation class
  //=================================================================

private:

  // Decides wether a radiative transition is possible and, if it is,
  // returns the identity of the starting shell for the transition
  G4int SelectTypeOfTransition(G4int Z, G4int shellId);
  
  // Generates a particle from a radiative transition and returns it
  G4DynamicParticle* GenerateFluorescence(G4int Z, G4int shellId, 
					  G4int provShellId);
 
  // Generates a particle from a non-radiative transition and returns it
  G4DynamicParticle* GenerateAuger(G4int Z, G4int shellId);

  //SI
  //Auger cascade by Burkhant Suerfu on March 24 2015 (Bugzilla 1727)
  //Generates auger electron cascade.
  G4DynamicParticle* GenerateAuger(G4int Z, G4int shellId, G4int& newAugerShellId);
  //ENDSI
  
  // copy constructor and hide assignment operator
  G4UAtomicDeexcitation(G4UAtomicDeexcitation &);
  G4UAtomicDeexcitation & operator=(const G4UAtomicDeexcitation &right);

  G4AtomicTransitionManager* transitionManager;
 
  // Data member which stores the shells to be filled by 
  // the radiative transition
  G4int newShellId;

  G4double minGammaEnergy;
  G4double minElectronEnergy;

  // Data member wich stores the id of the shell where is the vacancy 
  // left from the Auger electron
  G4int augerVacancyId;

  // Data member for the calculation of the proton and alpha ionisation XS

  G4VhShellCrossSection* PIXEshellCS;
  G4VhShellCrossSection* anaPIXEshellCS;
  G4VhShellCrossSection* ePIXEshellCS;
  G4EmCorrections*       emcorr;

  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;

  //SI
  //Auger cascade by Burkhant Suerfu on March 24 2015 (Bugzilla 1727)
  //Data member to keep track of cascading vacancies.
  std::vector<int> vacancyArray;
  //ENDSI
};

#endif




