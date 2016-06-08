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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VEvaporationFactory.cc,v 1.1 2002/06/06 18:04:27 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4VEvaporationFactory.hh"

const G4VEvaporationFactory & 
G4VEvaporationFactory::operator=(const G4VEvaporationFactory & val)
{
  G4Exception("G4VEvaporationFactory::operator= meant to not be accessable.");
  return *this;
}

G4bool 
G4VEvaporationFactory::operator==(const G4VEvaporationFactory & val) const
{
  G4Exception("G4VEvaporationFactory::operator== meant to not be accessable.");
  return false;
}

G4bool 
G4VEvaporationFactory::operator!=(const G4VEvaporationFactory & val) const
{
  G4Exception("G4VEvaporationFactory::operator!= meant to not be accessable.");
  return true;
}




G4VEvaporationFactory::~G4VEvaporationFactory()
{
  if (_channel != 0)
    G4std::for_each(_channel->begin(), _channel->end(), 
		    DeleteFragment());
  delete _channel;
}


G4std::vector<G4VEvaporationChannel*> * 
G4VEvaporationFactory::GetChannel()
{
  // Lazy initialization
  if (_channel == 0)
    _channel = CreateChannel();
  return _channel;
}



