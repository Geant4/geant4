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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: FissionLDP.hh,v 1.1 2003-10-08 12:32:17 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4FissionLevelDensityParameter_h_
#define G4FissionLevelDensityParameter_h_ 1


#include "G4VLevelDensityParameter.hh"
#include "MyLDP.hh"


class FissionLDP : public G4VLevelDensityParameter
{
public:
  FissionLDP() {};
  virtual ~FissionLDP() {};

private:  
  FissionLDP(const FissionLDP &right);

  const FissionLDP & operator=(const FissionLDP &right);
  G4bool operator==(const FissionLDP &right) const;
  G4bool operator!=(const FissionLDP &right) const;
  
public:
  G4double LevelDensityParameter(const G4int A,const G4int Z,const G4double U) const;

  
private:
  
  MyLDP theEvaporationLevelDensityParameter;

};


#endif
