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
#ifndef G4BatemanParameters_h
#define G4BatemanParameters_h 1

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4BatemanParameters.hh                                            //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   15 December 2015                                                  //
//  Description: renamed and streamlined version of F. Lei and P. Truscott's  //
//               class G4RadioactiveDecayRate.  This class contains the decay //
//               times and coefficients of the extended Bateman equations for //
//               the progeny of a given nuclide (Z, A, E).                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "globals.hh"
#include <vector>

class G4BatemanParameters
{

public:
  G4BatemanParameters();

  virtual ~G4BatemanParameters();

  // copy constructor and assignment operator
  G4BatemanParameters(const G4BatemanParameters &);
  G4BatemanParameters& operator=(const G4BatemanParameters&);

  // equality operators
  G4int operator==(const G4BatemanParameters& right) const
    {return (this == &right);};
  G4int operator!=(const G4BatemanParameters& right) const
    {return (this != &right);};

  void DumpInfo();

  inline G4int GetZ() const {return Z;}
  inline G4int GetA() const {return A;}
  inline G4double GetE() const {return E;}
  inline G4int GetGeneration() const {return generation;}
  inline std::vector<G4double> GetAcoefficients() const
    {return Acoeffs;}
  inline std::vector<G4double> GetTaus() const {return taus;}

  inline void SetZ(G4int value) {Z = value;}
  inline void SetA(G4int value) {A = value;}
  inline void SetE(G4double value) {E = value;}
  inline void SetGeneration(G4int value) {generation = value;}
  inline void SetAcoefficients(std::vector<G4double> value)
    {Acoeffs = value;}
  inline void SetTaus(std::vector<G4double> value) {taus = value;}

  void SetParameters(G4int /*Z*/, G4int /*A*/, G4double /*E*/, G4int /*G*/,
                     std::vector<G4double> /*Coeffs*/,
                     std::vector<G4double> /*taus*/);

private:
  G4int Z;
  G4int A;
  G4double E;
  G4int generation;
  std::vector<G4double> Acoeffs;
  std::vector<G4double> taus;

};
#endif
