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
#include "exrdmEnergyDeposition.hh"


//
// Default constructor
//
exrdmEnergyDeposition::exrdmEnergyDeposition()
{;}
//
// Specific constructor
//
exrdmEnergyDeposition::exrdmEnergyDeposition( G4double energy,
				    G4double time,
                                    G4double weight )
  : Energy(energy),
    Time(time),
    Weight(weight)
{;}


//
// Copy constructor
//
exrdmEnergyDeposition::exrdmEnergyDeposition( const exrdmEnergyDeposition &right )
  : Energy(right.Energy),
    Time(right.Time),
    Weight(right.Weight)
{;}

//
// Destructor
//
exrdmEnergyDeposition::~exrdmEnergyDeposition() {;}

//
// Equivalence operator
//
G4bool exrdmEnergyDeposition::operator==( const exrdmEnergyDeposition &right ) const
{
  return Time == right.Time;
}

//
// Order operators
//
G4bool exrdmEnergyDeposition::operator<( const exrdmEnergyDeposition &right ) const
{
  return Time < right.Time;
}

G4bool exrdmEnergyDeposition::operator<=( const exrdmEnergyDeposition &right ) const
{
  return Time <= right.Time;
}

