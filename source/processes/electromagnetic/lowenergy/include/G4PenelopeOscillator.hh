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
  G4PenelopeOscillator();
  G4PenelopeOscillator(const G4PenelopeOscillator&);
  
  ~G4PenelopeOscillator(){;};

  //I need to overload the following operators: > < == =
  G4PenelopeOscillator& operator=(const G4PenelopeOscillator&);
  int operator==(const G4PenelopeOscillator&) const;
  int operator>(const G4PenelopeOscillator&) const;
  int operator<(const G4PenelopeOscillator&) const;
  
  //Setters and getters
  G4double GetHartreeFactor() {return hartreeFactor;};
  void SetHartreeFactor(G4double hf) {hartreeFactor = hf;};
  //
  G4double GetIonisationEnergy() {return ionisationEnergy;};
  void SetIonisationEnergy(G4double ie) {ionisationEnergy = ie;};
  //
  G4double GetResonanceEnergy() const {return resonanceEnergy;};
  void SetResonanceEnergy(G4double re) {resonanceEnergy = re;};
  //
  G4double GetOscillatorStrength() {return oscillatorStrength;};
  void SetOscillatorStrength(G4double ostr) {oscillatorStrength=ostr;};
  //
  G4int GetShellFlag() {return shellFlag;};
  void SetShellFlag(G4int theflag) {shellFlag=theflag;};
  //
  G4double GetParentZ() {return parentZ;};
  void SetParentZ(G4double parZ) {parentZ = parZ;};
  //
  G4int GetParentShellID() {return parentShellID;};
  void SetParentShellID(G4int psID) {parentShellID = psID;};
  //
  G4double GetCutoffRecoilResonantEnergy(){return cutoffRecoilResonantEnergy;};
  void SetCutoffRecoilResonantEnergy(G4double ene){cutoffRecoilResonantEnergy = ene;};

private:
  G4double hartreeFactor;
  G4double ionisationEnergy;
  G4double resonanceEnergy;
  G4double oscillatorStrength;
  G4int shellFlag;
  G4double parentZ;
  G4int parentShellID; //necessary for interface to AtomicDeexcitationManager
  G4double cutoffRecoilResonantEnergy;
};

struct G4PenelopeOscillatorResEnergyComparator
{
public:
  int operator()(const G4PenelopeOscillator& left,
		 const G4PenelopeOscillator& right) 
  {return ((left.GetResonanceEnergy() < right.GetResonanceEnergy()) ? 1 : 0);};
};

#endif

