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
// GEANT4 physics class: G4HadronInelasticDataSet -- header file
// F.W. Jones, TRIUMF, 19-MAY-98
//
// Class Description
// Baseline data-set for  hadron inelastic cross-section. This does not need to
// be registered, but provides part of the general cross-section baseline 
// Class Description - End

#ifndef G4HadronInelasticDataSet_h
#define G4HadronInelasticDataSet_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4HadronCrossSections.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"


class G4HadronInelasticDataSet : public G4VCrossSectionDataSet
{
  public:

    G4HadronInelasticDataSet(); 

    virtual ~G4HadronInelasticDataSet()
    {}


    G4bool IsApplicable(const G4DynamicParticle* aParticle,
                        const G4Element* anElement)
    {
      return theHadronCrossSections->IsApplicable(aParticle, anElement);
    } 


    G4bool IsIsoApplicable(const G4DynamicParticle* aParticle,
                           G4int ZZ, G4int AA)
    {
      return theHadronCrossSections->IsApplicable(aParticle, ZZ, AA);
    }

    G4double GetCrossSection(const G4DynamicParticle* aParticle,
                             const G4Element* anElement, 
                             G4double /*aTemperature*/)
    {
      return theHadronCrossSections->GetInelasticCrossSection(aParticle,
                                                              anElement);
    }

    G4double GetZandACrossSection(const G4DynamicParticle* aParticle,
                                  G4int ZZ, G4int AA,
                                  G4double /*aTemperature*/)
    {
      return theHadronCrossSections->GetInelasticCrossSection(aParticle, ZZ, AA);
    }

    void BuildPhysicsTable(const G4ParticleDefinition&)
    {}

    void DumpPhysicsTable(const G4ParticleDefinition&);

  private:

    G4HadronCrossSections* theHadronCrossSections;
};

#endif
