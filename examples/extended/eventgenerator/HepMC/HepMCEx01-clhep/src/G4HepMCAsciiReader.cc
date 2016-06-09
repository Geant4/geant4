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
// ====================================================================
//
//   G4HepMCAsciiReader.cc
//   $Id: G4HepMCAsciiReader.cc,v 1.2 2006/06/29 17:09:16 gunter Exp $
//
// ====================================================================

#ifndef WIN32         // Temporarly disabled on Windows, until CLHEP
                      // will support the HepMC module

#include "G4HepMCAsciiReader.hh"
#include "G4HepMCAsciiReaderMessenger.hh"

//#include "g4std/iostream"

////////////////////////////////////////
G4HepMCAsciiReader::G4HepMCAsciiReader()
  :  filename("xxx.dat"), verbose(0)
////////////////////////////////////////
{
  messenger= new G4HepMCAsciiReaderMessenger(this);
}

/////////////////////////////////////////
G4HepMCAsciiReader::~G4HepMCAsciiReader()
/////////////////////////////////////////
{
  delete messenger;
}

/////////////////////////////////////
void G4HepMCAsciiReader::Initialize()
/////////////////////////////////////
{
  from.close();
  from.open(filename.c_str(), std::ios::in);
}

/////////////////////////////////////////////////////////
HepMC::GenEvent* G4HepMCAsciiReader::GenerateHepMCEvent()
/////////////////////////////////////////////////////////
{
  HepMC::GenEvent* evt= HepMC::readGenEvent(from);
  if(evt==0) {
    return 0; // no more event
    from.close();
  }

  if(verbose>0) evt-> print();
    
  return evt;
}

#endif
