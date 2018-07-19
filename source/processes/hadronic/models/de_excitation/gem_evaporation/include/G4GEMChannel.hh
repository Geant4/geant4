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
// $Id: G4GEMChannel.hh 107060 2017-11-01 15:00:04Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
#ifndef G4GEMChannel_h
#define G4GEMChannel_h 1

#include "G4VEvaporationChannel.hh"
#include "G4GEMProbability.hh"
#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4NucleiProperties.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"

class G4Pow;
class G4PairingCorrection;
class G4VCoulombBarrier;

class G4GEMChannel : public G4VEvaporationChannel
{
public:

  explicit G4GEMChannel(G4int theA, G4int theZ, const G4String & aName,
			G4GEMProbability * aEmissionStrategy);

  // destructor
  virtual ~G4GEMChannel();
    
  virtual G4double GetEmissionProbability(G4Fragment* theNucleus);

  virtual G4Fragment* EmittedFragment(G4Fragment* theNucleus);

  virtual void Dump() const;

  inline void SetLevelDensityParameter(G4VLevelDensityParameter * aLevelDensity)
  {
    if (MyOwnLevelDensity) { delete theLevelDensityPtr; }
    theLevelDensityPtr = aLevelDensity;
    MyOwnLevelDensity = false;
  }

private: 

  // Samples fragment kinetic energy.
  G4double SampleKineticEnergy(const G4Fragment & fragment);

  G4GEMChannel(const G4GEMChannel & right) = delete;  
  const G4GEMChannel & operator=(const G4GEMChannel & right) = delete;
  G4bool operator==(const G4GEMChannel & right) const = delete;
  G4bool operator!=(const G4GEMChannel & right) const = delete;
  
  // Data Members ************
  // This data member define the channel. 
  // They are intializated at object creation (constructor) time.
    
  // Atomic Number
  G4int A;
    
  // Charge
  G4int Z;

  G4double EvaporatedMass;
  G4double ResidualMass;

  G4Pow* fG4pow;
    
  // For evaporation probability calcualtion
  G4GEMProbability * theEvaporationProbabilityPtr;
    
  // For Level Density calculation
  G4bool MyOwnLevelDensity;
  G4VLevelDensityParameter * theLevelDensityPtr;
    
  // For Coulomb Barrier calculation
  G4VCoulombBarrier * theCoulombBarrierPtr;
  G4double CoulombBarrier;

  G4PairingCorrection* pairingCorrection;
    
  //---------------------------------------------------
    
  // These values depend on the nucleus that is being evaporated.
  // They are calculated through the Initialize method which takes as parameters 
  // the atomic number, charge and excitation energy of nucleus.
    
  // Residual Atomic Number
  G4int ResidualA;

  // Residual Charge
  G4int ResidualZ;
    
  // Emission Probability
  G4double EmissionProbability;

  // Maximal Kinetic Energy that can be carried by fragment
  G4double MaximalKineticEnergy;
};


#endif
