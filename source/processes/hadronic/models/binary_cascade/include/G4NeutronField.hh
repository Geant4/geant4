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
//      File name:     G4NeutronField.hh
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#ifndef G4NeutronField_h
#define  G4NeutronField_h 1

#include "G4VNuclearField.hh"
#include "G4V3DNucleus.hh"
#include "G4FermiMomentum.hh"
#include "G4VNuclearDensity.hh"
#include <vector>

class G4NeutronField: public G4VNuclearField
{
public:
  G4NeutronField(G4V3DNucleus * nucleus);
  virtual ~G4NeutronField();

private:

  G4NeutronField(const  G4NeutronField &right);
  const G4NeutronField & operator=(const G4NeutronField & right);
  int operator==(const G4NeutronField & right) const;
  int operator!=(const G4NeutronField & right) const;

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

  G4FermiMomentum theFermi;
  G4double theA;
  G4double theZ;
  const G4VNuclearDensity * theDensity;
  G4double theR;
  
  std::vector<G4double> theFermiMomBuffer;
};

#endif


