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
// ********************************************************************//
//
// $Id: G4coutDestination.hh 66241 2014-05-16 00:00:42Z adotti $
//
//
#include "G4coutIdDestination.hh"
#include "G4strstreambuf.hh"
#include "G4AutoLock.hh"

G4coutIdDestination::G4coutIdDestination(const G4int& idx,
                                          std::ostream& cout,
                                          std::ostream& cerr)
: G4coutDestination(),
    finalcout(cout),
    finalcerr(cerr),
    id(idx),
    useBuffer(false)
{
    G4coutbuf.SetDestination(this);
    G4cerrbuf.SetDestination(this);
}

G4coutIdDestination::~G4coutIdDestination()
{
    if (useBuffer) DumpBuffer();
}

void G4coutIdDestination::EnableBuffering(G4bool flag)
{
    useBuffer = flag;
}
namespace  {
    G4Mutex coutm = G4MUTEX_INITIALIZER;
}

void G4coutIdDestination::DumpBuffer()
{
    G4AutoLock l(&coutm);
    finalcout<<"====================="<<std::endl;
    finalcout<<"cout buffer for worker with ID:"<<id<<std::endl;
    finalcout<<cout_buffer.str()<<std::endl;
    finalcerr<<"====================="<<std::endl;
    finalcerr<<"cerr buffer for worker with ID:"<<id<<std::endl;
    finalcerr<<cerr_buffer.str()<<std::endl;
    finalcerr<<"====================="<<std::endl;
}

G4int G4coutIdDestination::ReceiveG4cout(const G4String& msg)
{
    if ( useBuffer )
        cout_buffer<<msg;
    else
        finalcout<<"G4Worker"<<id<<" > "<<msg;
    return 0;
}
G4int G4coutIdDestination::ReceiveG4cerr(const G4String& msg)
{
    if ( useBuffer )
        cerr_buffer<<msg;
    else
        finalcerr<<"G4Worker"<<id<<" > "<<msg;
    return 0;
}
