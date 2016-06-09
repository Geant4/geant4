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
#include "G4HETCEmissionFactory.hh"

#include "G4HETCNeutron.hh"
#include "G4HETCProton.hh"
#include "G4HETCDeuteron.hh"
#include "G4HETCTriton.hh"
#include "G4HETCHe3.hh"
#include "G4HETCAlpha.hh"


const G4HETCEmissionFactory & G4HETCEmissionFactory::
operator=(const G4HETCEmissionFactory & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4HETCEmissionFactory::operator= meant to not be accessable.");
  return *this;
}

G4bool G4HETCEmissionFactory::
operator==(const G4HETCEmissionFactory & ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4HETCEmissionFactory::operator== meant to not be accessable.");
  return false;
}

G4bool G4HETCEmissionFactory::
operator!=(const G4HETCEmissionFactory & ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4HETCEmissionFactory::operator!= meant to not be accessable.");
  return true;
}


std::vector<G4VPreCompoundFragment*> *  G4HETCEmissionFactory::
CreateFragmentVector()
{
  std::vector<G4VPreCompoundFragment*> * theFragVector = 
    new std::vector<G4VPreCompoundFragment*>;
  theFragVector->reserve(6);

  // neutron
  theFragVector->push_back(new G4HETCNeutron());
  // proton
  theFragVector->push_back(new G4HETCProton());
  // deuterium
  theFragVector->push_back(new G4HETCDeuteron());
  // triton
  theFragVector->push_back(new G4HETCTriton());
  // helium3
  theFragVector->push_back(new G4HETCHe3());
  // alpha
  theFragVector->push_back(new G4HETCAlpha());

  return theFragVector;
}
