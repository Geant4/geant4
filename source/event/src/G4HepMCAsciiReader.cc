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

// ====================================================================
//
//   G4HepMCAsciiReader.cc
//   $Id: G4HepMCAsciiReader.cc,v 1.2 2002-12-03 14:44:32 gcosmo Exp $
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
  from.open(filename.c_str(), G4std::ios::in);
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
