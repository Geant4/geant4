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
// GEANT4 physics class: G4HadronInelasticDataSet -- header file
// F.W. Jones, TRIUMF, 19-MAY-98
//
// Modified: V.Ivanchenko
//
// Class Description:
// Baseline data-set for  hadron inelastic cross-section. This does not need to
// be registered, but provides part of the general cross-section baseline 

#ifndef G4HadronInelasticDataSet_h
#define G4HadronInelasticDataSet_h 1

#include "G4VCrossSectionDataSet.hh"

class G4ParticleDefinition;
class G4NistManager;
class G4HadronCrossSections;

class G4HadronInelasticDataSet : public G4VCrossSectionDataSet
{
public:

  G4HadronInelasticDataSet(const G4String& name = "GheishaInelastic"); 

  virtual ~G4HadronInelasticDataSet();

  virtual void CrossSectionDescription(std::ostream&) const;

  virtual G4bool
  IsElementApplicable(const G4DynamicParticle* aParticle, G4int /*Z*/,
                      const G4Material*);

  virtual G4double
  GetElementCrossSection(const G4DynamicParticle* aParticle, G4int Z, 
			 const G4Material*);

private:

  G4int theZ;
  G4double fInelasticXS;
  G4double fKinEnergy;
  const G4ParticleDefinition* fParticle;
  G4HadronCrossSections* fGheishaXS;
  G4NistManager* fNIST;
};

#endif
