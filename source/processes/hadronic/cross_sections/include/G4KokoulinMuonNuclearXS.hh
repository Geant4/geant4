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
// Modified:
//
// 19 Aug 2011, V.Ivanchenko move to new design and make x-section per element
//

//
// Description: use Kokoulin's parameterized calculation of virtual 
//              photon production cross section and conversion to
//              real photons.

#ifndef G4KokoulinMuonNuclearXS_h
#define G4KokoulinMuonNuclearXS_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4PhysicsVector.hh"

const G4int MAXZMUN = 93;

class G4PhysicsVector;

class G4KokoulinMuonNuclearXS : public G4VCrossSectionDataSet
{
public:

  G4KokoulinMuonNuclearXS();
  virtual ~G4KokoulinMuonNuclearXS();

  static const char* Default_Name() {return "KokoulinMuonNuclearXS";}

  virtual void CrossSectionDescription(std::ostream&) const;

  G4bool IsElementApplicable(const G4DynamicParticle* particle, 
			     G4int Z, const G4Material*);

  G4double GetElementCrossSection(const G4DynamicParticle* particle,
				  G4int Z, const G4Material*);

  void BuildPhysicsTable(const G4ParticleDefinition&);

  void BuildCrossSectionTable();

  G4double
  ComputeDDMicroscopicCrossSection(G4double incidentKE, G4double Z,
                                   G4double A, G4double epsilon);

private:

  G4double
  ComputeMicroscopicCrossSection(G4double incidentKE, G4double A);

  G4KokoulinMuonNuclearXS & operator=(const G4KokoulinMuonNuclearXS &right);
  G4KokoulinMuonNuclearXS(const G4KokoulinMuonNuclearXS&);

  static G4PhysicsVector* theCrossSection[MAXZMUN];

  G4double LowestKineticEnergy;
  G4double HighestKineticEnergy;
  G4int    TotBin;
  G4double CutFixed;
  G4bool   isInitialized;
  G4bool   isMaster;
};

#endif
