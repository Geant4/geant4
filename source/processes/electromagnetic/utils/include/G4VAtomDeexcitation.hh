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
// $Id: G4VAtomDeexcitation.hh,v 1.1.2.1 2010/04/06 09:05:17 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-03-patch-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VAtomDeexcitation
//
// Author:        Alfonso Mantero & Vladimir Ivanchenko
//
// Creation date: 30.06.2009
//
// Modifications:
//
// Class Description:
//
// Abstract interface to energy loss models

// -------------------------------------------------------------------
//

#ifndef G4VAtomDeexcitation_h
#define G4VAtomDeexcitation_h 1

#include "globals.hh"
#include <vector>

class G4AtomicShell;
class G4ParticleDefinition;
class G4DynamicParticle;

class G4VAtomDeexcitation {

  G4VAtomDeexcitation(const G4String& pname = "");

  virtual ~G4VAtomDeexcitation();

  //========== initialization ==========

  virtual void PreparePhysicsTable(const G4ParticleDefinition&);
  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  // PIXE model name
  inline void SetPIXECrossSectionModel(const G4String&);
  inline const G4String& PIXECrossSectionModel() const;

  // Activation of deexcitation per detector region
  void SetFluorescenceActiveRegion(const G4Region* region = 0);
  void SetAugerActiveRegion(const G4Region* region = 0);
  void SetPIXECrossSectionActiveRegion(const G4Region* region = 0); 

  void SetFluorescenceActiveRegion(const G4String& rname = "");
  void SetAugerActiveRegion(const G4String& rname = "");
  void SetPIXECrossSectionActiveRegion(const G4String& rname = ""); 

  //========== Run time methods ==========

  // Check if deexcitation is active for a given geometry volume
  G4bool CheckFluorescenceActiveRegion(G4int coupleIndex);

  // Check if deexcitation is active for a given geometry volume
  G4bool CheckPIXEActiveRegion(G4int coupleIndex);

  // Get atomic shell by shell index, used by discrete processes 
  // (for example, photoelectric), when shell vacancy sampled by the model
  const G4AtomicShell* GetAtomicShell(G4int Z, G4int ShellIndex);

  // selection of random shell for ionisation process
  virtual const G4AtomicShell* SelectRandomShell(const G4DynamicParticle*, 
						 G4int Z);

  // generation of deexcitation for given atom and shell vacancy
  virtual void GenerateParticles(std::vector<G4DynamicParticle*>*, 
				 const G4AtomicShell*, G4int Z) = 0;

  // access or compute PIXE cross section 
  virtual G4double GetPIXECrossSection (const G4ParticleDefinition*, 
					G4int Z, G4double kinE) = 0;

  // calculate PIXE cross section from the models
  virtual G4double CalculatePIXECrossSection(const G4ParticleDefinition*,
					     G4int Z, G4double kinE) = 0;

  // Sampling of PIXE for ionisation processes
  virtual void 
  AlongStepDeexcitation(std::vector<G4DynamicParticle*>* secVect, 
			const G4DynamicParticle* icidentParticle, 
			const G4MaterialCutsCouple*, 
			G4double trueStepLenght, 
			G4double eLoss) = 0;



private:

  // copy constructor and hide assignment operator
  G4VAtomDeexcitation(G4VAtomDeexcitation &);
  G4VAtomDeexcitation & operator=(const G4VAtomDeexcitation &right);

  G4String namePIXE;

};

inline 
void G4VAtomDeexcitation::SetPIXECrossSectionModel(const G4String& n)
{
  namePIXE = n;
}

inline 
const G4String& G4VAtomDeexcitation::PIXECrossSectionModel() const
{
  return namePIXE;
}

#endif

