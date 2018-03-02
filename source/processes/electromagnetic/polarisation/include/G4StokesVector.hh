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
// $Id: G4StokesVector.hh 108502 2018-02-15 15:41:45Z gcosmo $
//
// GEANT4 Class header file
//
// File name:     G4StokesVector
//
// Author:        Andreas Schaelicke
//
// Creation date: 01.05.2005
//
// Modifications:
// 27-07-06 added some test routines (P.Starovoitov)
// 25-08-06 modified name of test routines (A.Schaelicke)
//
// Class Description:
//
// Provides Stokesvector representation employed in implementation of
// polarized processes.
//
// aim:
//   - store three components of a stokesvector
//   - distinguish between boson or fermion state (different transformations)
//   - provide unique definition of reference frame (cf. G4PolarizationHelper)
//

#ifndef G4StokesVector_h
#define G4StokesVector_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"


class G4StokesVector: public G4ThreeVector
{
 public:
  // standard vectors:
  static const G4StokesVector ZERO;
  static const G4StokesVector P1;
  static const G4StokesVector P2;
  static const G4StokesVector P3;
  static const G4StokesVector M1;
  static const G4StokesVector M2;
  static const G4StokesVector M3;
public:
  G4StokesVector();
  G4StokesVector(const G4ThreeVector & v);
  ~G4StokesVector() = default;

  G4bool IsZero() const; 

  inline G4double p1() const { return x(); }
  inline G4double p2() const { return y(); }
  inline G4double p3() const { return z(); }

  inline G4double Transverse() const { return perp(); } 

  inline G4ThreeVector PolSqr() const { 
    return G4ThreeVector(x()*x(),y()*y(),z()*z()); 
  }
  inline G4ThreeVector PolSqrt() const { 
    return G4ThreeVector(std::sqrt(x()),std::sqrt(y()),std::sqrt(z())); 
  }
  G4ThreeVector PolError(const G4StokesVector & sum2, long n);

  // Ratio of 3-vectors.
  G4ThreeVector PolDiv( const G4StokesVector & );

  inline void SetPhoton() { isPhoton=true; }

  void RotateAz(G4ThreeVector nInteractionFrame, 
		G4ThreeVector particleDirection);
  void InvRotateAz(G4ThreeVector nInteractionFrame, 
		   G4ThreeVector particleDirection);
  void RotateAz(G4double cosphi, G4double sinphi);
  G4double GetBeta();

  void DiceUniform();
  void DiceP1();
  void DiceP2();
  void DiceP3();

  void FlipP3();
private:
  G4bool isPhoton;
};


#endif
