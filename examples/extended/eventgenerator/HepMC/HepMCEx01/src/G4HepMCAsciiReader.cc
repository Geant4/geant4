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
//   $Id: G4HepMCAsciiReader.cc,v 1.1 2002-04-29 20:44:12 asaim Exp $
//
// ====================================================================
#include "G4HepMCAsciiReader.hh"
#include "G4HepMCAsciiReaderMessenger.hh"

#include "g4std/iostream"
#include "g4std/fstream"

////////////////////////////////////////
G4HepMCAsciiReader::G4HepMCAsciiReader()
  :  verbose(0), filename("xxx.dat")
////////////////////////////////////////
{
  asciiInput= new HepMC::IO_Ascii(filename.c_str(), ios::in);

  messenger= new G4HepMCAsciiReaderMessenger(this);
}

/////////////////////////////////////////
G4HepMCAsciiReader::~G4HepMCAsciiReader()
/////////////////////////////////////////
{
  delete asciiInput;
  delete messenger;
}

/////////////////////////////////////
void G4HepMCAsciiReader::Initialize()
/////////////////////////////////////
{
  delete asciiInput;

  asciiInput= new HepMC::IO_Ascii(filename.c_str(), ios::in);
  asciiInput-> print();
}

/////////////////////////////////////////////////////////
HepMC::GenEvent* G4HepMCAsciiReader::GenerateHepMCEvent()
/////////////////////////////////////////////////////////
{
  HepMC::GenEvent* evt= asciiInput-> read_next_event();
  if(!evt) return 0; // no more event

  if(verbose>0) evt-> print();
    
  return evt;
}

