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
// $Id: G4LPhysicsFreeVector.hh 98864 2016-08-15 11:53:26Z gcosmo $
//
// 
// ------------------------------------------------------------------
//
// Class G4LPhysicsFreeVector -- header file
//
// Class description:
//
// Derived from base class G4PhysicsVector
// This is a free vector for Low Energy Physics cross section data.
// The class name includes an "L" to distinguish it from other groups
// who may wish to implement a free vector in a different way.
// A subdivision method is used to find the energy|momentum bin.

// F.W. Jones, TRIUMF, 04-JUN-96
// 11-NOV-00 H.Kurashige: use STL vector for dataVector and binVector
// 02-APR-08 A.Bagulya: use GetValue() from base class
// 02-OCT-13  V.Ivanchenko : Remove FindBinLocation method
//
// ------------------------------------------------------------------

#ifndef G4LPhysicsFreeVector_h
#define G4LPhysicsFreeVector_h 1

#include "G4PhysicsFreeVector.hh"

class G4LPhysicsFreeVector : public G4PhysicsFreeVector  
{

public: // with description

  G4LPhysicsFreeVector();
  // the vector will be filled from external file using Retrieve method

  G4LPhysicsFreeVector(size_t length, G4double emin=0.0, G4double emax=0.0);
  // the vector with 'length' elements will be filled using PutValues method 
  // by default the vector is initialized with zeros

  virtual ~G4LPhysicsFreeVector();

  inline void PutValues(size_t index, G4double e, G4double dataValue);
  // user code is responsible for correct filling of all elements
};

inline
void G4LPhysicsFreeVector::PutValues(size_t index, G4double e, G4double value)
{
  G4PhysicsFreeVector::PutValue(index, e, value);
}

#endif
