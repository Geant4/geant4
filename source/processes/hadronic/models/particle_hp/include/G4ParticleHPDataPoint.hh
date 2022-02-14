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
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPDataPoint_h
#define G4ParticleHPDataPoint_h 1

#include "globals.hh"

class G4ParticleHPDataPoint
{
  public:
  
  G4ParticleHPDataPoint(){energy = 0; xSec = 0;}
  G4ParticleHPDataPoint(G4double e, G4double x){ energy = e; xSec = x;}
  
  G4ParticleHPDataPoint & operator= (const G4ParticleHPDataPoint & aSet) = default;

//  ~G4ParticleHPDataPoint(){}
  
  inline G4double GetEnergy() const   {return energy;}
  inline G4double GetXsection() const {return xSec;}
  
  inline void SetEnergy(G4double e)  {energy = e;}
  inline void SetXsection(G4double x){xSec = x;}
  
  inline G4double GetX() const {return energy;}
  inline G4double GetY() const {return xSec;}
  
  inline void SetX(G4double e)  {energy = e;}
  inline void SetY(G4double x)  {xSec = x;}
  
  inline void SetData(G4double e, G4double x) {energy = e; xSec = x;}
  
  private:
  
  G4double energy;
  G4double xSec;
};

#endif
