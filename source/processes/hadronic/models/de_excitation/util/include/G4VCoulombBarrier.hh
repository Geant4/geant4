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
// $Id: G4VCoulombBarrier.hh,v 1.2 2005/06/04 13:29:20 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)

#ifndef G4VCoulombBarrier_h
#define G4VCoulombBarrier_h 1

#include "globals.hh"


class G4VCoulombBarrier
{
public:
  G4VCoulombBarrier(const G4int anA, const G4int aZ);
  virtual ~G4VCoulombBarrier() {};

protected:
  G4VCoulombBarrier() : theA(1),theZ(0) {};

private:
  G4VCoulombBarrier(const G4VCoulombBarrier & right);

  const G4VCoulombBarrier & operator=(const G4VCoulombBarrier & right);
  G4bool operator==(const G4VCoulombBarrier & right) const;
  G4bool operator!=(const G4VCoulombBarrier & right) const;
  
public:
  virtual G4double GetCoulombBarrier(const G4int ARes, const G4int ZRes, 
				     const G4double U) const = 0;
					
  G4int GetA(void) const {return theA;}
  G4int GetZ(void) const {return theZ;}
					
private: 
	
  G4int theA;
  G4int theZ;

};

#endif
