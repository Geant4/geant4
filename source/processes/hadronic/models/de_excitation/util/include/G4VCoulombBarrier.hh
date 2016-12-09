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
// $Id: G4VCoulombBarrier.hh 100690 2016-10-31 11:25:43Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)

#ifndef G4VCoulombBarrier_h
#define G4VCoulombBarrier_h 1

#include "globals.hh"

class G4VCoulombBarrier
{
public:

  explicit G4VCoulombBarrier(G4int anA, G4int aZ);
  virtual ~G4VCoulombBarrier();

  virtual G4double GetCoulombBarrier(G4int ARes, G4int ZRes, 
				     G4double U) const = 0;

  virtual G4double BarrierPenetrationFactor(G4int Eexc) const = 0;

  void SetParameters(G4double rho, G4double r0); 
					
  inline G4int GetA(void) const { return theA; }
  inline G4int GetZ(void) const { return theZ; }
  inline G4double GetRho(void) const { return theRho; }
  inline G4double GetR0(void)  const { return theR0; }

private:

  G4VCoulombBarrier(const G4VCoulombBarrier & right) = delete;
  const G4VCoulombBarrier & operator=(const G4VCoulombBarrier & right) = delete;
  G4bool operator==(const G4VCoulombBarrier & right) const = delete;
  G4bool operator!=(const G4VCoulombBarrier & right) const = delete;
	
  G4int theA;
  G4int theZ;
  G4double theRho;
  G4double theR0;
};

#endif
