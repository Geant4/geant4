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

// by V. Lara

#include "G4VPreCompoundEmissionFactory.hh"

const G4VPreCompoundEmissionFactory & 
G4VPreCompoundEmissionFactory::operator=(const G4VPreCompoundEmissionFactory & )
{
  G4Exception("G4VPreCompoundEmissionFactory::operator= meant to not be accessable.");
  return *this;
}

G4bool 
G4VPreCompoundEmissionFactory::operator==(const G4VPreCompoundEmissionFactory & ) const
{
  G4Exception("G4VPreCompoundEmissionFactory::operator== meant to not be accessable.");
  return false;
}

G4bool 
G4VPreCompoundEmissionFactory::operator!=(const G4VPreCompoundEmissionFactory & ) const
{
  G4Exception("G4VPreCompoundEmissionFactory::operator!= meant to not be accessable.");
  return true;
}




G4VPreCompoundEmissionFactory::~G4VPreCompoundEmissionFactory()
{
  if (_fragvector != 0)
    std::for_each(_fragvector->begin(), _fragvector->end(), 
		    DeleteFragment());
  delete _fragvector;
}


std::vector<G4VPreCompoundFragment*> * 
G4VPreCompoundEmissionFactory::GetFragmentVector()
{
  // Lazy initialization
  if (_fragvector == 0)
    _fragvector = CreateFragmentVector();
  return _fragvector;
}

