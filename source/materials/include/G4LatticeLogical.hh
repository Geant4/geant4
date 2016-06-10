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
/// \file materials/include/G4LatticeLogical.hh
/// \brief Definition of the G4LatticeLogical class
//
// $Id: G4LatticeLogical.hh 76883 2013-11-18 12:50:08Z gcosmo $
//
// 20131114  Add verbosity for diagnostic output
// 20131115  Expose maximum array dimensions for use by LatticeReader

#ifndef G4LatticeLogical_h
#define G4LatticeLogical_h

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>


class G4LatticeLogical {
public:
  G4LatticeLogical();
  virtual ~G4LatticeLogical();

  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  G4bool LoadMap(G4int, G4int, G4int, G4String);
  G4bool Load_NMap(G4int, G4int, G4int, G4String);

  // Get group velocity magnitude for input polarization and wavevector
  virtual G4double MapKtoV(G4int, const G4ThreeVector& ) const;

  // Get group velocity direction (unit vector) for input polarization and K
  virtual G4ThreeVector MapKtoVDir(G4int, const G4ThreeVector& ) const;

  // Dump structure in format compatible with reading back
  void Dump(std::ostream& os) const;
  void DumpMap(std::ostream& os, G4int pol, const G4String& name) const;
  void Dump_NMap(std::ostream& os, G4int pol, const G4String& name) const;

public:
  void SetDynamicalConstants(G4double Beta, G4double Gamma,
			     G4double Lambda, G4double Mu) {
    fBeta=Beta; fGamma=Gamma; fLambda=Lambda; fMu=Mu;
  }

  void SetScatteringConstant(G4double b) { fB=b; }
  void SetAnhDecConstant(G4double a) { fA=a; }
  void SetLDOS(G4double LDOS) { fLDOS=LDOS; }
  void SetSTDOS(G4double STDOS) { fSTDOS=STDOS; }
  void SetFTDOS(G4double FTDOS) { fFTDOS=FTDOS; }

  G4double GetBeta() const { return fBeta; }
  G4double GetGamma() const { return fGamma; }
  G4double GetLambda() const { return fLambda; }
  G4double GetMu() const { return fMu; }
  G4double GetScatteringConstant() const { return fB; }
  G4double GetAnhDecConstant() const { return fA; }
  G4double GetLDOS() const { return fLDOS; }
  G4double GetSTDOS() const { return fSTDOS; }
  G4double GetFTDOS() const { return fFTDOS; }

public:
  enum { MAXRES=322 };			    // Maximum map resolution (bins)
  
private:
  G4int verboseLevel;			    // Enable diagnostic output

  G4double fMap[3][MAXRES][MAXRES];	    // map for group velocity scalars
  G4ThreeVector fN_map[3][MAXRES][MAXRES];  // map for direction vectors

  G4int fVresTheta; //velocity  map theta resolution (inclination)
  G4int fVresPhi;   //velocity  map phi resolution  (azimuth)
  G4int fDresTheta; //direction map theta resn
  G4int fDresPhi;   //direction map phi resn 

  G4double fA;       //Scaling constant for Anh.Dec. mean free path
  G4double fB;       //Scaling constant for Iso.Scat. mean free path
  G4double fLDOS;    //Density of states for L-phonons
  G4double fSTDOS;   //Density of states for ST-phonons
  G4double fFTDOS;   //Density of states for FT-phonons
  G4double fBeta, fGamma, fLambda, fMu; //dynamical constants for material
};

// Write lattice structure to output stream

inline std::ostream& 
operator<<(std::ostream& os, const G4LatticeLogical& lattice) {
  lattice.Dump(os);
  return os;
}

#endif
