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
// $Id: G4ITManager_iterator.cc 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITManager.hh"

// Description of G4ITManager<IT>::iterator

G4VITManager::iterator::iterator(G4ITBox* box)
{
    if(box)
    {
        fCurrentBox = box ;
    }
    else 
    {
        fCurrentBox = 0;	
    }

    fNextIT = 0;
}
//______________________________________________________________________

G4bool G4VITManager::iterator::begin()
{
    if(fCurrentBox)
    {
        fNextIT = fCurrentBox->GetFirstIT();

        if(fNextIT)
        {
//            G4cout << "G4VITManager::allbox_iterator, fNextIT :" <<  fNextIT -> GetName()<< G4endl;
            return true;
        }
    }

    return false;
}
//______________________________________________________________________

G4VITManager::iterator& G4VITManager::iterator::operator++(G4int)
{
    if(fCurrentBox)
    {
        if(fNextIT)
        {
            fNextIT = fNextIT->GetNext();
        }
    }
    return *this;
}

//______________________________________________________________________

G4bool G4VITManager::iterator::end()
{
    if(fNextIT) return false;
    return true ;
}

//______________________________________________________________________

G4VITManager::iterator& G4VITManager::iterator::operator= (const iterator& i)
{
    if(this != &i)
    {
        fNextIT = i.fNextIT;
        fCurrentBox = i.fCurrentBox;
    }
    return *this;
}

//______________________________________________________________________

G4IT* G4VITManager::iterator::operator*()
{
    return fNextIT;
}

//______________________________________________________________________

G4ITBox* G4VITManager::iterator::GetBox()
{
    return fCurrentBox;
}

//______________________________________________________________________

void G4VITManager::iterator::PrintNext() const
{
    if(fNextIT)
    {
        if(fNextIT->GetTrack())
        {
                G4cout
                    << fNextIT->GetTrack()->GetTrackID()
                    << "\t"
                    << fNextIT->GetName()
                    << G4endl;
        }
    }
    else
    {
        G4cout<<"fNextIT = 0"<<G4endl;
    }
}

