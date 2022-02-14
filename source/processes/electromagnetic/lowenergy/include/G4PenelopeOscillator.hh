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
// Author: Luciano Pandola
//
// History:
// -----------
// 18 Dec 2008   L. Pandola   First implementation
//                            
//
// -------------------------------------------------------------------
//
// Class description:
// Class designed to contain data of atomic oscillators that are used by 
// several Penelope models. Objects of G4PenelopeOscillators are managed 
// by a dedicated G4PenelopeOscillatorManager.
// -------------------------------------------------------------------

#ifndef G4PENELOPEOSCILLATOR_HH
#define G4PENELOPEOSCILLATOR_HH 1

#include "globals.hh"

class G4PenelopeOscillator 
{
 public:
  explicit G4PenelopeOscillator();
  G4PenelopeOscillator(const G4PenelopeOscillator&);
  
  ~G4PenelopeOscillator(){;};

  //I need to overload the following operators: > < == =
  G4PenelopeOscillator& operator=(const G4PenelopeOscillator&);
  G4bool operator==(const G4PenelopeOscillator&) const;
  G4bool operator>(const G4PenelopeOscillator&) const;
  G4bool operator<(const G4PenelopeOscillator&) const;
  
  //Setters and getters
  G4double GetHartreeFactor() {return fHartreeFactor;};
  void SetHartreeFactor(G4double hf) {fHartreeFactor = hf;};
  //
  G4double GetIonisationEnergy() {return fIonisationEnergy;};
  void SetIonisationEnergy(G4double ie) {fIonisationEnergy = ie;};
  //
  G4double GetResonanceEnergy() const {return fResonanceEnergy;};
  void SetResonanceEnergy(G4double re) {fResonanceEnergy = re;};
  //
  G4double GetOscillatorStrength() {return fOscillatorStrength;};
  void SetOscillatorStrength(G4double ostr) {fOscillatorStrength=ostr;};
  //
  G4int GetShellFlag() {return fShellFlag;};
  void SetShellFlag(G4int theflag) {fShellFlag=theflag;};
  //
  G4double GetParentZ() {return fParentZ;};
  void SetParentZ(G4double parZ) {fParentZ = parZ;};
  //
  G4int GetParentShellID() {return fParentShellID;};
  void SetParentShellID(G4int psID) {fParentShellID = psID;};
  //
  G4double GetCutoffRecoilResonantEnergy(){return fCutoffRecoilResonantEnergy;};
  void SetCutoffRecoilResonantEnergy(G4double ene){fCutoffRecoilResonantEnergy = ene;};

private:
  G4double fHartreeFactor;
  G4double fIonisationEnergy;
  G4double fResonanceEnergy;
  G4double fOscillatorStrength;
  G4double fParentZ;  
  G4double fCutoffRecoilResonantEnergy;
  G4int fParentShellID; //necessary for interface to AtomicDeexcitationManager
  G4int fShellFlag;
};

struct G4PenelopeOscillatorResEnergyComparator
{
public:
  G4int operator()(const G4PenelopeOscillator& left,
		 const G4PenelopeOscillator& right) 
  {return ((left.GetResonanceEnergy() < right.GetResonanceEnergy()) ? 1 : 0);};
};

#endif

