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

  G4double GetEnergy() {return Energy;};
  G4double GetTime() {return Time;};
  G4double GetWeight() {return Weight;};
  // Accessors

  private:

    G4double Energy;  
    G4double Time;    
    G4double Weight;
};
#endif



