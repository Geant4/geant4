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
// $Id: G4VNuclearDensity.hh,v 1.3.8.1 2001/06/28 19:13:29 gunter Exp $
// GEANT4 tag $Name:  $
//
#ifndef G4VNuclearDensity_h
#define G4VNuclearDensity_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"


class G4VNuclearDensity 
{

  public:
    G4VNuclearDensity();
    virtual ~G4VNuclearDensity();
    
    inline G4double GetDensity(G4ThreeVector aPosition)
    {
	return rho0*GetRelativeDensity(aPosition);
    };
    
    virtual G4double GetRelativeDensity(G4ThreeVector aPosition) = 0;
    virtual G4double GetRadius(const G4double maxRelativeDenisty) = 0;
    virtual G4double GetDeriv(const G4ThreeVector & point) = 0;    

  protected:    
    inline void Setrho0(G4double arho0) { rho0=arho0; };
    inline G4double Getrho0() { return rho0; };
   
  private:
  
    G4double rho0;
};

#endif

