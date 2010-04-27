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
//
// $Id: G4VEvaporationFactory.cc,v 1.7 2010-04-27 11:43:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4VEvaporationFactory.hh"

G4VEvaporationFactory::G4VEvaporationFactory() : _channel(0) 
{}

G4VEvaporationFactory::~G4VEvaporationFactory()
{
  if (_channel != 0)
    std::for_each(_channel->begin(), _channel->end(), 
		    DeleteFragment());
  delete _channel;
}


std::vector<G4VEvaporationChannel*> * 
G4VEvaporationFactory::GetChannel()
{
  // Lazy initialization
  if (_channel == 0)
    _channel = CreateChannel();
  return _channel;
}



