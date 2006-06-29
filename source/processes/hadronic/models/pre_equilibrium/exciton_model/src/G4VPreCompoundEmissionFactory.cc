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

// by V. Lara

#include "G4VPreCompoundEmissionFactory.hh"
#include "G4HadronicException.hh"

const G4VPreCompoundEmissionFactory & 
G4VPreCompoundEmissionFactory::operator=(const G4VPreCompoundEmissionFactory & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4VPreCompoundEmissionFactory::operator= meant to not be accessable.");
  return *this;
}

G4bool 
G4VPreCompoundEmissionFactory::operator==(const G4VPreCompoundEmissionFactory & ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4VPreCompoundEmissionFactory::operator== meant to not be accessable.");
  return false;
}

G4bool 
G4VPreCompoundEmissionFactory::operator!=(const G4VPreCompoundEmissionFactory & ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4VPreCompoundEmissionFactory::operator!= meant to not be accessable.");
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

