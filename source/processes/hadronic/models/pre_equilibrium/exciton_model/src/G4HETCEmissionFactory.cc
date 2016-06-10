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
// $Id: G4HETCEmissionFactory.cc 68028 2013-03-13 13:48:15Z gcosmo $
//
// by V. Lara
//
// Modified:
// 23.08.2010 V.Ivanchenko general cleanup, move constructor and destructor 
//            the source
//

#include "G4HETCEmissionFactory.hh"

#include "G4HETCNeutron.hh"
#include "G4HETCProton.hh"
#include "G4HETCDeuteron.hh"
#include "G4HETCTriton.hh"
#include "G4HETCHe3.hh"
#include "G4HETCAlpha.hh"

G4HETCEmissionFactory::G4HETCEmissionFactory()
{}

G4HETCEmissionFactory::~G4HETCEmissionFactory()
{}

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
  // alpha
  theFragVector->push_back(new G4HETCAlpha());
  // triton
  theFragVector->push_back(new G4HETCTriton());
  // helium3
  theFragVector->push_back(new G4HETCHe3());

  return theFragVector;
}
