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
// $Id$
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITManager.hh"

using namespace std;

__thread auto_ptr<G4AllITManager> *G4AllITManager::fInstance_G4MT_TLS_ = 0;

G4AllITManager::G4AllITManager()
{ if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ;
    fVerbose = 0 ;
}

G4AllITManager* G4AllITManager::Instance()
{  ;;;   if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ; auto_ptr<G4AllITManager> &fInstance = *fInstance_G4MT_TLS_;  ;;;  
    if(fInstance.get() == 0) fInstance = auto_ptr<G4AllITManager>(new G4AllITManager());
    return fInstance.get() ;
}

void G4AllITManager::DeleteInstance()
{  ;;;   if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ; auto_ptr<G4AllITManager> &fInstance = *fInstance_G4MT_TLS_;  ;;;  
    fInstance.reset();
}

G4AllITManager::~G4AllITManager()
{  ;;;   if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ; auto_ptr<G4AllITManager> &fInstance = *fInstance_G4MT_TLS_;  ;;;  
    std::map<G4ITType, G4VITManager*>::iterator it ;
    std::map<G4ITType, G4VITManager*>::iterator it_tmp ;

    for(it = fITSubManager.begin(); it!=fITSubManager.end() ; )
    {
        if(it->second) delete it->second;
        it_tmp = it;
        it++;
        fITSubManager.erase(it_tmp);
    }
    fInstance.release();
}

void  G4AllITManager::UpdatePositionMap()
{ if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ;
    std::map<G4ITType, G4VITManager*>::iterator it = fITSubManager.begin() ;

    for(; it!=fITSubManager.end() ; it++)
    {
        it->second->UpdatePositionMap();
    }
}

void  G4AllITManager::CreateTree()
{ if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ;
    std::map<G4ITType, G4VITManager*>::iterator it = fITSubManager.begin() ;

    for(; it!=fITSubManager.end() ; it++)
    {
        it->second->CreateTree();
    }
}

template<typename T>  G4ITManager<T>* G4AllITManager::Instance()
{  ;;;   if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ; auto_ptr<G4AllITManager> &fInstance = *fInstance_G4MT_TLS_;  ;;;  
    return G4ITManager<T>::Instance();
}

G4VITManager* G4AllITManager::GetInstance(G4ITType type)
{ if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ;
    map<G4ITType, G4VITManager*>::iterator it = fITSubManager.find(type);

    if(it == fITSubManager.end()) return 0;

    return it->second;
}

void G4AllITManager::RegisterManager(G4VITManager* manager)
{ if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ;
    fITSubManager[manager->GetITType()] = manager;
}

G4ITBox* G4AllITManager::GetBox(const G4Track* track)
{ if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ;
    map<G4ITType, G4VITManager*>::iterator it = fITSubManager.find(GetIT(track)->GetITType());

    if(it == fITSubManager.end()) return 0;

    return it->second->GetBox(track);
}

void G4AllITManager::Push(G4Track* track)
{ if (!fInstance_G4MT_TLS_) fInstance_G4MT_TLS_ = new auto_ptr<G4AllITManager> (0) ;
    fITSubManager[GetIT(track)->GetITType()]->Push(track);
}
