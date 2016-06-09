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
// $Id: G4ITManager_allbox_iterator.cc 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITManager.hh"

// Definitions of G4VITManager::allbox_iterator

G4VITManager::allbox_iterator::allbox_iterator(G4ITType type) : iterator(0)
{
    fInstance = G4AllITManager::Instance()->GetInstance(type);
    fCurrentBox = fInstance->GetFirstBox();
    fNextIT = 0 ;
}

G4VITManager::allbox_iterator::allbox_iterator(G4VITManager* man) : iterator(0)
{
    fInstance = man;
    fCurrentBox = fInstance->GetFirstBox();
    fNextIT = 0 ;
}

//______________________________________________________________________

G4VITManager::allbox_iterator& G4VITManager::allbox_iterator::operator++(G4int)
{	
    if(fNextIT)
    {
        fNextIT = fNextIT -> GetNext();
        if(fNextIT == 0)
        {
            while (fCurrentBox && fCurrentBox != fInstance->GetLastBox() && fNextIT == 0)
            {
                fCurrentBox = fCurrentBox->GetNextBox();
                if(fCurrentBox)
                {
                    fNextIT = fCurrentBox->GetFirstIT();
//                    if(fNextIT)
//                        G4cout << "allbox_iterator::operator++ -- : "
//                               << fNextIT -> GetName()
//                               << " (" << fNextIT->GetTrack()->GetTrackID() << ") "
//                               << G4endl;
                }
            }
        }
//        else
//        {
//            G4cout << "allbox_iterator::operator++ -- : "
//                   << fNextIT -> GetName()
//                   << " (" << fNextIT->GetTrack()->GetTrackID() << ") "
//                   << G4endl;
//        }
    }
    return *this;
}

//______________________________________________________________________
G4VITManager::allbox_iterator& G4VITManager::allbox_iterator::operator= (const allbox_iterator& i)
{
    if(this != &i)
    {
        fNextIT = i.fNextIT;
        fCurrentBox = i.fCurrentBox;
    }
    return *this;
}

