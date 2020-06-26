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
// G4VRangeToEnergyConverter
//
// Class description:
//
// Base class for Range to Energy Converters.
// Cut in energy corresponding to given cut value in range
// is calculated for a material by using Convert() method.

// Author: H.Kurashige, 05 October 2002 - First implementation
// --------------------------------------------------------------------
#ifndef G4VRangeToEnergyConverter_hh
#define G4VRangeToEnergyConverter_hh 1

#include <cmath>
#include <vector>

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"

class G4PhysicsLogVector;

class G4VRangeToEnergyConverter
{
  public:

    G4VRangeToEnergyConverter();
      // Constructor

    G4VRangeToEnergyConverter(const G4VRangeToEnergyConverter& r);
      // Copy constructor

    G4VRangeToEnergyConverter& operator=(const G4VRangeToEnergyConverter &r);
      // Assignment operator

    virtual ~G4VRangeToEnergyConverter();
      // Destructor

    G4bool operator==(const G4VRangeToEnergyConverter& r) const;
    G4bool operator!=(const G4VRangeToEnergyConverter& r) const;
      // Equality operators

    virtual G4double Convert(G4double rangeCut, const G4Material* material);
      // Calculate energy cut from given range cut for the material

    static void SetEnergyRange(G4double lowedge, G4double highedge);
      // Set energy range for all particle type

    static G4double GetLowEdgeEnergy();
    static G4double GetHighEdgeEnergy();
      // Get energy range for all particle type

    static G4double GetMaxEnergyCut();
    static void SetMaxEnergyCut(G4double value);
      // Get/set max cut energy for all particle type
  
    inline const G4ParticleDefinition* GetParticleType() const;
      // Return pointer to the particle type which this converter takes care of

    const  G4PhysicsTable* GetLossTable() const;   
      // theLossTable is a collection of loss vectors for all elements.
      // Each loss vector has energy loss values (cross-section values
      // for neutral particles) which are calculated by
      // ComputeLoss(G4double AtomicNumber, G4double KineticEnergy).
      // ComputeLoss method is pure virtual and should be provided
      // for each particle type

    virtual void Reset();
      // Reset Loss Table and Range Vectors

    inline void SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;
      // control flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

  protected:

    virtual void BuildLossTable();

    virtual G4double ComputeLoss(G4double AtomicNumber,
                                 G4double KineticEnergy) = 0;

    // ------------- Range Table --------------------------------------

    using G4LossVector = G4PhysicsLogVector;
    using G4RangeVector = G4PhysicsLogVector;
    using G4LossTable = G4PhysicsTable;

    virtual void BuildRangeVector(const G4Material* aMaterial,
                                  G4RangeVector* rangeVector);

    G4double ConvertCutToKineticEnergy(G4RangeVector* theRangeVector,
                                       G4double       theCutInLength, 
                                       std::size_t    materialIndex ) const;
  protected:

    static G4double LowestEnergy, HighestEnergy;
    static G4double MaxEnergyCut; 
    G4double fMaxEnergyCut = 0.0;
   
    const G4ParticleDefinition* theParticle = nullptr;
    G4LossTable* theLossTable = nullptr;
    G4int NumberOfElements = 0;
  
    const G4int TotBin = 300;

    std::vector< G4RangeVector* > fRangeVectorStore;   

  private:

    G4int verboseLevel = 1;
};

// ------------------
// Inline methods
// ------------------

inline 
void G4VRangeToEnergyConverter::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline 
G4int G4VRangeToEnergyConverter::GetVerboseLevel() const
{
  return verboseLevel;
}

inline 
const G4ParticleDefinition* G4VRangeToEnergyConverter::GetParticleType() const
{
  return theParticle;
}

#endif
