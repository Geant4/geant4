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
// V. Ivanchenko 20 October 2023 
//
// Cross-section data set for a high precision (based on evaluated data
// libraries) description of light ion inelastic interactions below 200 MeV.
//

#ifndef G4ParticleHPInelasticXS_h
#define G4ParticleHPInelasticXS_h 1

#include "G4CrossSectionHP.hh"
#include "G4ParticleDefinition.hh"
#include <fstream>

class G4ParticleHPInelasticXS final : public G4CrossSectionHP
{
  public:
    explicit G4ParticleHPInelasticXS(const G4ParticleDefinition*);

    ~G4ParticleHPInelasticXS() override = default;

    void CrossSectionDescription(std::ostream&) const final;

    G4ParticleHPInelasticXS & operator=(const G4ParticleHPInelasticXS &right) = delete;
    G4ParticleHPInelasticXS(const G4ParticleHPInelasticXS&) = delete;

 private:

    const G4ParticleDefinition* part;
};

#endif
