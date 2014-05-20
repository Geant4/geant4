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
//
#ifndef G4Nucleon_h
#define G4Nucleon_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4Nucleon ----------------
//             by Gunter Folger, May 1998.
//       class for a nucleon (inside a 3D Nucleus)
// ------------------------------------------------------------

#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4AntiProton.hh"    // Uzhi Feb. 2011
#include "G4AntiNeutron.hh"   // Uzhi Feb. 2011 

#include "G4VKineticNucleon.hh"

//#include "G4VSplitableHadron.hh"
class G4VSplitableHadron;

class G4Nucleon : public G4VKineticNucleon
{

  public:
      G4Nucleon();
      ~G4Nucleon();

      inline int operator==(const G4Nucleon &right) const;
      inline int operator!=(const G4Nucleon &right) const;
      G4Nucleon& operator=(const G4Nucleon& right);

  public:

    inline void SetPosition(G4ThreeVector & aPosition) {thePosition = aPosition;}
    virtual inline const G4ThreeVector & GetPosition()  const {return thePosition;}

    inline void SetMomentum(G4LorentzVector & aMomentum) {theMomentum = aMomentum;}
    inline const G4LorentzVector& GetMomentum()  const {return theMomentum;}
    virtual inline const G4LorentzVector & Get4Momentum()  const {return theMomentum;}

    inline void SetBindingEnergy(G4double anEnergy) {theBindingE = anEnergy;}
    inline G4double GetBindingEnergy()  const {return theBindingE;}

    inline void SetParticleType(G4Proton * aProton) {theParticleType = aProton;}
    inline void SetParticleType(G4Neutron *aNeutron){theParticleType = aNeutron;}

    inline void SetParticleType(G4AntiProton * aAntiProton) {theParticleType =aAntiProton;} //VU
    inline void SetParticleType(G4AntiNeutron *aAntiNeutron){theParticleType =aAntiNeutron;}//VU


    inline  const G4ParticleDefinition* GetParticleType() const {return theParticleType;}
    virtual const G4ParticleDefinition* GetDefinition() const {return theParticleType;}
    
    inline void Boost(const G4ThreeVector & beta){ theMomentum.boost(beta); } 
           void Boost(const G4LorentzVector & aMomentum);

    inline void Hit(G4VSplitableHadron *aHit) { theSplitableHadron=aHit;}
//    inline void Hit(G4int ) { isHit=true;}    
    inline void Hit(G4int ) 
    { 
      theSplitableHadron=reinterpret_cast<G4VSplitableHadron *>(1111); 
    }
    inline G4VSplitableHadron * GetSplitableHadron() const { return theSplitableHadron;}
    inline G4bool AreYouHit() const {  return theSplitableHadron!=0;}

  private:

    G4ThreeVector thePosition;
    G4LorentzVector theMomentum;
    G4double theBindingE;
    const G4ParticleDefinition * theParticleType;
    G4VSplitableHadron * theSplitableHadron;


};

std::ostream & operator << (std::ostream &, const G4Nucleon&);

inline int G4Nucleon::operator==(const G4Nucleon &right) const
{
	return this==&right;
}
inline int G4Nucleon::operator!=(const G4Nucleon &right) const
{
	return this!=&right;
}

inline G4Nucleon& G4Nucleon::operator=(const G4Nucleon& right)
{
   if (this != &right)
   {
      thePosition=right.GetPosition();
      theMomentum=right.Get4Momentum();
      theBindingE=right.GetBindingEnergy();
      theParticleType=right.GetDefinition();
      theSplitableHadron=right.GetSplitableHadron();
   }
	return *this;
}

#endif


