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
#include "G4PreCompoundEmissionFactory.hh"

#include "G4PreCompoundNeutron.hh"
#include "G4PreCompoundProton.hh"
#include "G4PreCompoundDeuteron.hh"
#include "G4PreCompoundTriton.hh"
#include "G4PreCompoundHe3.hh"
#include "G4PreCompoundAlpha.hh"


const G4PreCompoundEmissionFactory & G4PreCompoundEmissionFactory::
operator=(const G4PreCompoundEmissionFactory & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4PreCompoundEmissionFactory::operator= meant to not be accessable.");
  return *this;
}

G4bool G4PreCompoundEmissionFactory::
operator==(const G4PreCompoundEmissionFactory & ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4PreCompoundEmissionFactory::operator== meant to not be accessable.");
  return false;
}

G4bool G4PreCompoundEmissionFactory::
operator!=(const G4PreCompoundEmissionFactory & ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4PreCompoundEmissionFactory::operator!= meant to not be accessable.");
  return true;
}


std::vector<G4VPreCompoundFragment*> *  G4PreCompoundEmissionFactory::
CreateFragmentVector()
{
  std::vector<G4VPreCompoundFragment*> * theFragVector = 
    new std::vector<G4VPreCompoundFragment*>();
  theFragVector->reserve(6);

  // neutron
  theFragVector->push_back(new G4PreCompoundNeutron());
  // proton
  theFragVector->push_back(new G4PreCompoundProton());
  // deuterium
  theFragVector->push_back(new G4PreCompoundDeuteron());
  // alpha
  theFragVector->push_back(new G4PreCompoundAlpha());
  // triton
  theFragVector->push_back(new G4PreCompoundTriton());
  // helium3
  theFragVector->push_back(new G4PreCompoundHe3());

  return theFragVector;
}
