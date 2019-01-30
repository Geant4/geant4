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
//      File name:     G4ProtonField.hh
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#ifndef G4ProtonField_h
#define  G4ProtonField_h 1

#include "G4VNuclearField.hh"
#include "G4V3DNucleus.hh"
#include "G4FermiMomentum.hh"
#include "G4VNuclearDensity.hh"
#include <vector>

class G4ProtonField: public G4VNuclearField
{
public:

  G4ProtonField(G4V3DNucleus * nucleus);
  virtual ~G4ProtonField();

private:
  G4ProtonField(const  G4ProtonField &right);
  const G4ProtonField & operator=(const G4ProtonField & right);
  G4bool operator==(const G4ProtonField & right) const;
  G4bool operator!=(const G4ProtonField & right) const;

public:

  virtual G4double GetField(const G4ThreeVector & aPosition);
  virtual G4double GetBarrier();

private:

  G4double GetDensity(const G4ThreeVector & aPosition)
  {
    return theDensity->GetDensity(aPosition);
  }
  G4double GetFermiMomentum(const G4double aDensity)
  {
    return theFermi.GetFermiMomentum(aDensity);
  }
  
  G4int theA;
  G4int theZ;
  G4double theBarrier;
  G4double theRadius;
  G4FermiMomentum theFermi;
  const G4VNuclearDensity * theDensity;

  std::vector<G4double> theFermiMomBuffer;
};

#endif


