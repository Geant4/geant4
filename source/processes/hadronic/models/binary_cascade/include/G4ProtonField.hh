//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
  int operator==(const G4ProtonField & right) const;
  int operator!=(const G4ProtonField & right) const;

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
  
  G4double theA;
  G4double theZ;
  G4double theBarrier;
  G4double theRadius;
  G4FermiMomentum theFermi;
  const G4VNuclearDensity * theDensity;

  std::vector<G4double> theFermiMomBuffer;
};

#endif


