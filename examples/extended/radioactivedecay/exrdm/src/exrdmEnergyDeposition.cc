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

