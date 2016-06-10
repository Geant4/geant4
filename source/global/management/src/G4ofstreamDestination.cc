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
// $Id: G4coutDestination.hh 66241 2014-05-16 00:00:42Z adotti $
//
//
// ---------------------------------------------------------------
// GEANT 4 class cc file
//
// G4ofstreamDestination.cc
//
// ---------------------------------------------------------------
#include "G4strstreambuf.hh"
#include "G4ofstreamDestination.hh"
#include <ios>

//TODO: make this that can be only cout / only cerr or both with a parameter,
//  default is both
G4ofstreamDestinationBase::G4ofstreamDestinationBase(const G4String& fn ,
                                             G4bool append)
    : G4coutDestination(),
    fileName(fn),
    appendFlag(append)
{
}

G4ofstreamDestinationBase::~G4ofstreamDestinationBase()
{
    Close();
}

void G4ofstreamDestinationBase::SetFileName(const G4String& fn , G4bool append )
{
    fileName = fn;
    appendFlag = append;
}

void G4ofstreamDestinationBase::Open()
{
    if ( ! g4file.is_open() )
    {
        std::ios::openmode mode = std::ios::out;
        if ( appendFlag )
            mode |= std::ios::app;
        g4file.open(fileName,mode);
    }
}

void G4ofstreamDestinationBase::Close()
{
    if ( g4file.is_open() )
        g4file.close();
}


//========================
//Concrete implementations
//========================
G4CoutToFile::G4CoutToFile(const G4String& defaultName , G4bool append) :
    G4ofstreamDestinationBase(defaultName,append)
{
    G4coutbuf.SetDestination(this);
}

G4int G4CoutToFile::ReceiveG4cout(const G4String& msg)
{
    if ( ! g4file.is_open() ) Open();
    g4file<<msg<<std::flush;
    return 0;
}

G4CerrToFile::G4CerrToFile(const G4String& defaultName, G4bool append) :
    G4ofstreamDestinationBase(defaultName,append)
{
    G4cerrbuf.SetDestination(this);
}

G4int G4CerrToFile::ReceiveG4cerr(const G4String& msg)
{
    if ( !g4file.is_open() ) Open();
    g4file<<msg<<std::flush;
    return 0;
}
