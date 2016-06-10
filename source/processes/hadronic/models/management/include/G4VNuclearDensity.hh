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
// $Id: G4VNuclearDensity.hh 66785 2013-01-12 15:10:13Z gcosmo $
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
    
    inline G4double GetDensity(const G4ThreeVector & aPosition) const
    {
	return rho0*GetRelativeDensity(aPosition);
    };
    
    virtual G4double GetRelativeDensity(const G4ThreeVector & aPosition) const = 0;
    virtual G4double GetRadius(const G4double maxRelativeDenisty) const = 0;
    virtual G4double GetDeriv(const G4ThreeVector & point) const = 0;    

  protected:    
    inline void Setrho0(G4double arho0) { rho0=arho0; };
    inline G4double Getrho0() const { return rho0; };
   
  private:
  
    G4double rho0;
};

#endif

