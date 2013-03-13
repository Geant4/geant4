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
#ifndef G4Fancy3DNucleus_h
#define G4Fancy3DNucleus_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4Fancy3DNucleus ----------------
//             by Gunter Folger, May 1998.
//       class for a 3D nucleus, arranging nucleons in space and momentum.
// ------------------------------------------------------------
// 20110805  M. Kelsey -- Remove C-style array (pointer) of G4Nucleons,
//		make vector a container of objects.  Move testSums,
//		places, momentum and fermiM to class data members for
//		reuse.  Remove args from ReduceSum(), use data members.

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4Nucleon.hh"		/* FIXME: This should be forward decl! */
#include "G4V3DNucleus.hh"
#include "G4VNuclearDensity.hh"
#include "G4FermiMomentum.hh"
#include <vector>

class G4Fancy3DNucleusHelper;

// to test if we can drop old interface for (A,Z), comment next line..
//#define NON_INTEGER_A_Z 1

class G4Fancy3DNucleus : public G4V3DNucleus
{

  public:
      G4Fancy3DNucleus();
      ~G4Fancy3DNucleus();

  private:
      G4Fancy3DNucleus(const G4Fancy3DNucleus &right);
      const G4Fancy3DNucleus & operator=(const G4Fancy3DNucleus &right);
      int operator==(const G4Fancy3DNucleus &right) const;
      int operator!=(const G4Fancy3DNucleus &right) const;
      

//  Implementation 
      void ChooseNucleons();
      void ChoosePositions();
      void ChooseFermiMomenta();
      G4double BindingEnergy();
      G4bool ReduceSum();

  public:
#if defined(NON_INTEGER_A_Z)
      void Init(G4double theA, G4double theZ);
#endif
      void Init(G4int theA, G4int theZ);
      G4bool StartLoop();
      G4Nucleon * GetNextNucleon();
      const std::vector<G4Nucleon> & GetNucleons();
      G4int GetMassNumber();
      G4double GetMass();
      G4int GetCharge();
      G4double GetNuclearRadius();
      G4double GetNuclearRadius(const G4double maxRelativeDensity);
      G4double GetOuterRadius();
      G4double AddExcitationEnergy(G4double);
      G4double GetExcitationEnergy();
      G4double CoulombBarrier();
      void DoLorentzBoost(const G4LorentzVector & theBoost);
      void DoLorentzBoost(const G4ThreeVector & theBeta);
      void DoLorentzContraction(const G4LorentzVector & theBoost);
      void DoLorentzContraction(const G4ThreeVector & theBeta);
      void CenterNucleons();
      void DoTranslation(const G4ThreeVector & theShift);
      const G4VNuclearDensity * GetNuclearDensity() const;
      void SortNucleonsIncZ();            // on increased Z-coordinates Uzhi 29.08.08
      void SortNucleonsDecZ();            // on decreased Z-coordinates Uzhi 29.08.08
      
  private:
  
  G4int myA;
  G4int myZ;
  std::vector<G4Nucleon> theNucleons;

  G4int currentNucleon;
  G4VNuclearDensity * theDensity;
  G4FermiMomentum theFermi;  
  const G4double nucleondistance;
  G4double excitationEnergy;
  
  std::vector<G4ThreeVector> places;		// For selecting locations
  std::vector<G4ThreeVector> momentum;		// For selecting nucleon motion
  std::vector<G4double> fermiM;
  std::vector<G4Fancy3DNucleusHelper> testSums;	// For sorting nucleon configs
};


inline G4int G4Fancy3DNucleus::GetCharge()
{
	return myZ;
}

inline G4int G4Fancy3DNucleus::GetMassNumber()
{
	return myA;
}
inline G4double G4Fancy3DNucleus::AddExcitationEnergy(G4double anE)
{
   excitationEnergy +=anE;
   return excitationEnergy;
}

inline G4double G4Fancy3DNucleus::GetExcitationEnergy()
{
   return excitationEnergy;
}

#endif
