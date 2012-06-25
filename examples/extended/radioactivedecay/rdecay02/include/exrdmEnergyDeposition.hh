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
/// \file radioactivedecay/rdecay02/include/exrdmEnergyDeposition.hh
/// \brief Definition of the exrdmEnergyDeposition class
//
#ifndef exrdmEnergyDeposition_h
#define exrdmEnergyDeposition_h 1

#include "globals.hh"

class exrdmEnergyDeposition
{
  public:   // with description

    exrdmEnergyDeposition();
    exrdmEnergyDeposition( const exrdmEnergyDeposition &right );
    exrdmEnergyDeposition( G4double, G4double, G4double );
    virtual ~exrdmEnergyDeposition();
         // Constructor and virtual destructor

    G4bool operator==(const exrdmEnergyDeposition &right) const ;
    G4bool operator< (const exrdmEnergyDeposition &right) const ;
    G4bool operator<=(const exrdmEnergyDeposition &right) const ;
  // Operators  

  G4double GetEnergy() {return fEnergy;};
  G4double GetTime() {return fTime;};
  G4double GetWeight() {return fWeight;};
  // Accessors

  private:

    G4double fEnergy;  
    G4double fTime;    
    G4double fWeight;
};
#endif



