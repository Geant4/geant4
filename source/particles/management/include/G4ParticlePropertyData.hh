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
// $Id: G4ParticlePropertyData.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: 
// ---------------- G4ParticlePropertyData ----------------
// first implementation by H Kurashige 9 June 2003
// Add   magnetic moment    by H Kurashige   Mar 2007
// ------------------------------------------------------------

#ifndef G4ParticlePropertyData_h
#define G4ParticlePropertyData_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4ParticlePropertyTable;
class G4ParticlePropertyData 
{
  // Class Description
  //  This class containes properties of a particle which are subset
  //  of properties in G4ParticleDefinition class.
  //  This class is used only for G4ParticlePropertyTable.
  //
  
  friend class  G4ParticlePropertyTable;

 public: // With Description
 
  G4ParticlePropertyData(const G4String& particleName = "");
  // The particle name should be assigned
  // This particle name can not be changed except for assignment operator

  G4ParticlePropertyData(const G4ParticlePropertyData &right);
 
  virtual ~G4ParticlePropertyData();
      
  G4ParticlePropertyData & operator=(const G4ParticlePropertyData &right);
  
  G4int operator==(const G4ParticlePropertyData &right) const;
  G4int operator!=(const G4ParticlePropertyData &right) const;

 public: // With Description
  // By these following Getxxxx methods, you can get values 
  // for members which can not be changed
  const G4String& GetParticleName() const { return theParticleName; }
  
  G4double GetPDGMass() const { return thePDGMass; }
  G4double GetPDGWidth() const { return thePDGWidth; } 
  G4double GetPDGCharge() const { return thePDGCharge; }
  
  G4int    GetPDGiSpin() const { return thePDGiSpin; }
  G4int    GetPDGiParity() const { return thePDGiParity; }
  G4int    GetPDGiConjugation() const { return thePDGiConjugation; }
  G4int    GetPDGiIsospin() const { return thePDGiIsospin; }
  G4int    GetPDGiIsospin3() const { return thePDGiIsospin3; }
  G4int    GetPDGiGParity() const { return thePDGiGParity; }
 
  G4double GetPDGMagneticMoment() const { return thePDGMagneticMoment; }
   
  G4int    GetLeptonNumber() const { return theLeptonNumber; }
  G4int    GetBaryonNumber() const { return theBaryonNumber; }
  
  G4int    GetPDGEncoding() const { return thePDGEncoding; }
  G4int    GetAntiPDGEncoding() const { return theAntiPDGEncoding; }
  
  G4int    GetQuarkContent(G4int flavor) const;
  G4int    GetAntiQuarkContent(G4int flavor) const;
  //  return the number of quark with flavor contained in this particle. 
  //  The value of flavor is assigned as follows 
  // 1:d, 2:u, 3:s, 4:c, 5:b, 6:t
  
  G4double GetPDGLifeTime() const { return thePDGLifeTime; }

  // SetXXX methods 
  void SetPDGMass(G4double newMass);
  void SetPDGWidth(G4double newWidth);
  void SetPDGCharge(G4double newCharge);
  
  void SetPDGiSpin(G4int newSpin);
  void SetPDGiParity(G4int newParity);
  void SetPDGiConjugation(G4int newConjugation);
  void SetPDGiIsospin(G4int newIsospin);
  void SetPDGiIsospin3(G4int newIsospin3);
  void SetPDGiGParity(G4int newGParity);

  void SetPDGMagneticMoment(G4double mageticMoment);
  
  void SetLeptonNumber(G4int newLeptonNumber);
  void SetBaryonNumber(G4int newBaryonNumber);
  
  void SetPDGEncoding(G4int newEncoding);
  void SetAntiPDGEncoding(G4int newAntiEncoding);
  
  void SetQuarkContent(G4int flavor, G4int newContent);
  void SetAntiQuarkContent(G4int flavor, G4int newContent);
  void SetPDGLifeTime(G4double newLifeTime); 

 public: // With Description
  void Print() const;
  //  Prints information of data members.

 public:
  void  SetVerboseLevel(G4int value);
  G4int GetVerboseLevel() const;
  // controle flag for output message
  //  0: Silent
  //  1: Warning message
  //  2: More


 private:
  G4String theParticleName;
  //  The name of the particle.
  
  G4double thePDGMass;
  //  The mass of the particle, in units of equivalent energy.
  
  G4double thePDGWidth;
  //  The decay width of the particle, usually the width of a
  //  Breit-Wigner function, assuming that you are near the
  //  mass center anyway. (in units of equivalent energy)
  
  G4double thePDGCharge;
  //  The charge of the particle.(in units of Coulomb)
  
  //   ---- following members are quantum number
  //         i.e. discrete numbers can be allowded
  //        So, you can defined only by using integer in constructor 

  G4int thePDGiSpin;
  //  The total spin of the particle, also often denoted as
  //  capital J, in units of 1/2.

  G4int thePDGiParity;
  //  The parity quantum number, in units of 1. If the parity
  //  is not defined for this particle, we will set this to 0.
  
  G4int thePDGiConjugation;
  //  This charge conjugation quantum number in units of 1.
  
  G4int thePDGiGParity;
  //  The value of the G-parity quantum number.
  
  G4int thePDGiIsospin;
  G4int thePDGiIsospin3;
  //  The isospin and its 3rd-component in units of 1/2.
  
  G4double thePDGMagneticMoment;
  //  The magnetic moment.

  G4int theLeptonNumber;
  //  The lepton quantum number.
  
  G4int theBaryonNumber;
  //  The baryon quantum number.
  
  G4int thePDGEncoding;
  //  The Particle Data Group integer identifier of this particle
  
  G4int theAntiPDGEncoding;
  //  The Particle Data Group integer identifier of the anti-particle
  
  G4double thePDGLifeTime;
  //  The Particle Life Time
  
  enum {NumberOfQuarkFlavor = 6};
  G4int  theQuarkContent[NumberOfQuarkFlavor];
  G4int  theAntiQuarkContent[NumberOfQuarkFlavor];
  //  the number of quark (minus Sign means anti-quark) contents

 private:
  G4bool fPDGMassModified;
  G4bool fPDGWidthModified;
  G4bool fPDGChargeModified;
  G4bool fPDGiSpinModified;
  G4bool fPDGiParityModified;
  G4bool fPDGiConjugationModified;
  G4bool fPDGiGParityModified;
  G4bool fPDGiIsospinModified;
  G4bool fPDGiIsospin3Modified;
  G4bool fPDGIsospinModified;
  G4bool fPDGIsospin3Modified;
  G4bool fPDGMagneticMomentModified;
  G4bool fLeptonNumberModified;
  G4bool fBaryonNumberModified;
  G4bool fPDGEncodingModified;
  G4bool fAntiPDGEncodingModified;
  G4bool fQuarkContentModified;
  G4bool fAntiQuarkContentModified;
  G4bool fPDGLifeTimeModified;

 private:
  G4int verboseLevel;

};

#include "G4ParticlePropertyData.icc"

#endif


