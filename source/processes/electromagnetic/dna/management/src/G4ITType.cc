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
// $Id: G4ITType.cc 80151 2014-04-03 09:42:22Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITType.hh"
#include "G4AutoLock.hh"

/*G4ThreadLocal*/ G4ITTypeManager* G4ITTypeManager::fgInstance = 0;
G4ThreadLocal G4ITTypeManager* G4ITTypeManager::fgInstance_local = 0;

G4Mutex deleteMutex=G4MUTEX_INITIALIZER;
G4Mutex ressourceMutex=G4MUTEX_INITIALIZER;

// static method
size_t G4ITType::size()
{
    return  G4ITTypeManager::Instance()->size();
}

G4ITType& G4ITType::operator=(const G4ITType & rhs)
{
    if (this == &rhs) return *this;
    fValue = rhs.fValue;
    return *this;
}

G4ITTypeManager*  G4ITTypeManager::Instance()
{
    if(fgInstance == 0)
    {
        fgInstance =  new G4ITTypeManager();
    }
    return fgInstance;
}

void G4ITTypeManager::DeleteInstance()
{
	G4AutoLock lock(&deleteMutex);
	if(fgInstance)
	{
		delete fgInstance ;
		fgInstance = 0;
	}
}

void G4ITTypeManager::ReserveRessource()
{
	G4AutoLock lock(&ressourceMutex);
	fRessource++;
}

void G4ITTypeManager::ReleaseRessource()
{
	G4AutoLock lock(&ressourceMutex);
	fRessource--;

	if(fRessource <= 0) DeleteInstance();
}

G4ITTypeManager::G4ITTypeManager()
{
	fLastType = 0;
	fRessource = 0;
}

G4ITTypeManager::~G4ITTypeManager()
{;}

size_t G4ITTypeManager::size() const
{
    return fLastType;
}

G4ITType G4ITTypeManager::NewType()
{
    G4ITType newType = fLastType;
    fLastType++;
    return newType;
}
