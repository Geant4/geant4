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
// G4ParticlePropertyData
//
// Class description:
//
// This class contains properties of a particle which are a subset
// of properties in G4ParticleDefinition.
// This class is used only for G4ParticlePropertyTable.

// Author: H.Kurashige, 9 June 2003
// --------------------------------------------------------------------
#ifndef G4ParticlePropertyData_hh
#define G4ParticlePropertyData_hh 1

#include "globals.hh"
#include "G4ios.hh"

class G4ParticlePropertyTable;

class G4ParticlePropertyData 
{
  friend class G4ParticlePropertyTable;

  public:
 
    G4ParticlePropertyData(const G4String& particleName = "");
      // The particle name should be assigned
      // The name cannot be changed except by assignment operator

    G4ParticlePropertyData(const G4ParticlePropertyData& right);
 
    virtual ~G4ParticlePropertyData();
      
    G4ParticlePropertyData& operator=(const G4ParticlePropertyData& right);

    G4bool operator==(const G4ParticlePropertyData& right) const;
    G4bool operator!=(const G4ParticlePropertyData& right) const;

    // With the following accessors, one can get values 
    // for members which cannot be changed

    const G4String& GetParticleName() const { return theParticleName; }
  
    G4double GetPDGMass() const { return thePDGMass; }
    G4double GetPDGWidth() const { return thePDGWidth; } 
    G4double GetPDGCharge() const { return thePDGCharge; }
  
    G4int GetPDGiSpin() const { return thePDGiSpin; }
    G4int GetPDGiParity() const { return thePDGiParity; }
    G4int GetPDGiConjugation() const { return thePDGiConjugation; }
    G4int GetPDGiIsospin() const { return thePDGiIsospin; }
    G4int GetPDGiIsospin3() const { return thePDGiIsospin3; }
    G4int GetPDGiGParity() const { return thePDGiGParity; }
 
    G4double GetPDGMagneticMoment() const { return thePDGMagneticMoment; }
   
    G4int GetLeptonNumber() const { return theLeptonNumber; }
    G4int GetBaryonNumber() const { return theBaryonNumber; }
  
    G4int GetPDGEncoding() const { return thePDGEncoding; }
    G4int GetAntiPDGEncoding() const { return theAntiPDGEncoding; }
  
    inline G4int GetQuarkContent(G4int flavor) const;
    inline G4int GetAntiQuarkContent(G4int flavor) const;
      // Return the number of quark with flavor contained in this particle. 
      // The value of flavor is assigned as follows: 
      // 1:d, 2:u, 3:s, 4:c, 5:b, 6:t
  
    G4double GetPDGLifeTime() const { return thePDGLifeTime; }

    // Modifiers

    inline void SetPDGMass(G4double newMass);
    inline void SetPDGWidth(G4double newWidth);
    inline void SetPDGCharge(G4double newCharge);
  
    inline void SetPDGiSpin(G4int newSpin);
    inline void SetPDGiParity(G4int newParity);
    inline void SetPDGiConjugation(G4int newConjugation);
    inline void SetPDGiIsospin(G4int newIsospin);
    inline void SetPDGiIsospin3(G4int newIsospin3);
    inline void SetPDGiGParity(G4int newGParity);

    inline void SetPDGMagneticMoment(G4double magneticMoment);
  
    inline void SetLeptonNumber(G4int newLeptonNumber);
    inline void SetBaryonNumber(G4int newBaryonNumber);
  
    inline void SetPDGEncoding(G4int newEncoding);
    inline void SetAntiPDGEncoding(G4int newAntiEncoding);
  
    inline void SetQuarkContent(G4int flavor, G4int newContent);
    inline void SetAntiQuarkContent(G4int flavor, G4int newContent);
    inline void SetPDGLifeTime(G4double newLifeTime); 

    void Print() const;
      // Prints information of data members

    inline void SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;
      // Control flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

  private:

    G4String theParticleName = "";
      // The name of the particle.
  
    G4double thePDGMass = 0.0;
      // The mass of the particle, in units of equivalent energy.
  
    G4double thePDGWidth = 0.0;
      // The decay width of the particle, usually the width of a
      // Breit-Wigner function, assuming that you are near the
      // mass center anyway (in units of equivalent energy).
  
    G4double thePDGCharge = 0.0;
      // The charge of the particle (in units of Coulomb).
  
    //   ---- following members are quantum number
    //        i.e. discrete numbers can be allowed
    //        So, you can define only by using integer in constructor 

    G4int thePDGiSpin = 0;
      // The total spin of the particle, also often denoted as
      // capital J, in units of 1/2.

    G4int thePDGiParity = 0;
      // The parity quantum number, in units of 1. If the parity
      // is not defined for this particle, we will set this to 0.
  
    G4int thePDGiConjugation = 0;
      // This charge conjugation quantum number in units of 1.
  
    G4int thePDGiGParity = 0;
      // The value of the G-parity quantum number.
  
    G4int thePDGiIsospin = 0;
    G4int thePDGiIsospin3 = 0;
      // The isospin and its 3rd-component in units of 1/2.
  
    G4double thePDGMagneticMoment = 0.0;
      // The magnetic moment.

    G4int theLeptonNumber = 0;
      // The lepton quantum number.
  
    G4int theBaryonNumber = 0;
      // The baryon quantum number.
  
    G4int thePDGEncoding = 0;
      // The Particle Data Group integer identifier of this particle
  
    G4int theAntiPDGEncoding = 0;
      // The Particle Data Group integer identifier of the anti-particle
  
    G4double thePDGLifeTime = -1.0;
      // The Particle Life Time
  
    enum { NumberOfQuarkFlavor = 6 };

    G4int theQuarkContent[NumberOfQuarkFlavor];
    G4int theAntiQuarkContent[NumberOfQuarkFlavor];
      // The number of quark (minus Sign means anti-quark) contents

  private:

    G4bool fPDGMassModified = false;
    G4bool fPDGWidthModified = false;
    G4bool fPDGChargeModified = false;
    G4bool fPDGiSpinModified = false;
    G4bool fPDGiParityModified = false;
    G4bool fPDGiConjugationModified = false;
    G4bool fPDGiGParityModified = false;
    G4bool fPDGiIsospinModified = false;
    G4bool fPDGiIsospin3Modified = false;
    G4bool fPDGIsospinModified = false;
    G4bool fPDGIsospin3Modified = false;
    G4bool fPDGMagneticMomentModified = false;
    G4bool fLeptonNumberModified = false;
    G4bool fBaryonNumberModified = false;
    G4bool fPDGEncodingModified = false;
    G4bool fAntiPDGEncodingModified = false;
    G4bool fQuarkContentModified = false;
    G4bool fAntiQuarkContentModified = false;
    G4bool fPDGLifeTimeModified = false;

    G4int verboseLevel = 1;
};

#include "G4ParticlePropertyData.icc"

#endif
