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
// $Id: G4PreCompoundEmissionFactory.cc 68028 2013-03-13 13:48:15Z gcosmo $
//

#include "G4PreCompoundEmissionFactory.hh"

#include "G4PreCompoundNeutron.hh"
#include "G4PreCompoundProton.hh"
#include "G4PreCompoundDeuteron.hh"
#include "G4PreCompoundTriton.hh"
#include "G4PreCompoundHe3.hh"
#include "G4PreCompoundAlpha.hh"

G4PreCompoundEmissionFactory::G4PreCompoundEmissionFactory()
{}

G4PreCompoundEmissionFactory::~G4PreCompoundEmissionFactory()
{}

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
