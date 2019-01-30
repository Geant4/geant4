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
// -------------------------------------------------------------------
//      GEANT 4 class header file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4KaonPlusField.hh
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#ifndef G4KaonPlusField_h
#define  G4KaonPlusField_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VNuclearField.hh"
#include "G4V3DNucleus.hh"

class G4KaonPlusField: public G4VNuclearField
{
public:
  G4KaonPlusField(G4V3DNucleus * nucleus, G4double coeff = 0.35*CLHEP::fermi);
  virtual ~G4KaonPlusField();

private:
  G4KaonPlusField(const  G4KaonPlusField &right);
  const G4KaonPlusField & operator=(const G4KaonPlusField & right);
  G4bool operator==(const G4KaonPlusField & right) const;
  G4bool operator!=(const G4KaonPlusField & right) const;

public:
  virtual G4double GetField(const G4ThreeVector & aPosition);
  virtual G4double GetBarrier();
  virtual G4double GetCoeff() { return theCoeff; }

private:
  G4double theCoeff;
};

#endif



