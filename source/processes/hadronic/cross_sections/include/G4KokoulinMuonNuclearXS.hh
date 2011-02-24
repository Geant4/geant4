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
// $Id: $
//
// Author:      D.H. Wright (SLAC)
// Date:        1 February 2011
//
// Description: use Kokoulin's parameterized calculation of virtual 
//              photon production cross section and conversion to
//              real photons.

#ifndef G4KokoulinMuonNuclearXS_h
#define G4KokoulinMuonNuclearXS_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
// #include "G4ParticleTable.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NucleiProperties.hh"
#include <vector>


class G4PhysicsTable;

class G4KokoulinMuonNuclearXS : public G4VCrossSectionDataSet
{
public:

  G4KokoulinMuonNuclearXS();
  virtual ~G4KokoulinMuonNuclearXS();


  G4bool IsApplicable(const G4DynamicParticle* particle, const G4Element* )
  {
    return IsIsoApplicable(particle, 0, 0);
  }

  G4bool IsIsoApplicable(const G4DynamicParticle* particle,
                         G4int /*ZZ*/, G4int /*AA*/)
  {
    G4bool result = false;
    if (particle->GetDefinition() == G4MuonMinus::MuonMinus() ||
        particle->GetDefinition() == G4MuonPlus::MuonPlus() ) result = true;

    return result;
  }


  G4double GetCrossSection(const G4DynamicParticle* particle, 
                           const G4Element* element, G4double temp = 0.);


  G4double GetZandACrossSection(const G4DynamicParticle* particle,
                                G4int ZZ, G4int AA, G4double /*aTemperature*/);

  void BuildCrossSectionTable();

  void BuildPhysicsTable(const G4ParticleDefinition&) {}

  void DumpPhysicsTable(const G4ParticleDefinition&) {}

  G4double
  ComputeDDMicroscopicCrossSection(G4double incidentKE, G4double,
                                   G4double AtomicWeight, G4double epsilon);

private:

  G4double
  ComputeMicroscopicCrossSection(G4double incidentKE,
                                 G4double AtomicNumber, G4double AtomicWeight);

  G4PhysicsTable* theCrossSectionTable;

  G4double LowestKineticEnergy;
  G4double HighestKineticEnergy;
  G4int TotBin;
  G4double CutFixed;

  std::map<G4int,G4Element*> zelMap;
};

#endif
