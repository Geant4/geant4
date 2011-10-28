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
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4AllITManager_h
#define G4AllITManager_h 1
#include <globals.hh>
#include <map>
#include <vector>
#include "G4ITType.hh"
#include "G4ThreeVector.hh"

class G4IT;
class G4VITManager;
class G4ITBox;
class G4Track;
template<typename T> class G4ITManager;
/*
Class Description :
*/
class G4AllITManager
{
private :
    static G4AllITManager* fInstance;
    G4AllITManager();

    std::map<G4ITType, G4VITManager*> fITSubManager ;

    int fVerbose ;

public :
    static G4AllITManager* Instance();
    ~G4AllITManager();

    template<typename T>  G4ITManager<T>* Instance();

    G4VITManager* GetInstance(G4ITType);
    G4ITBox* GetBox(const G4Track*);

    void RegisterManager(G4VITManager* manager);
    void Push(G4Track* track);

    void SetVerboseLevel(G4int level)
    {
        fVerbose = level;
    }
    G4int GetVerboseLevel()
    {
        return fVerbose;
    }

    void UpdatePositionMap();
    void CreateTree();

    template<typename T> std::vector<std::pair<G4IT*, double> >* FindNearest(const G4ThreeVector& pos, const T* it);
    template<typename T> std::vector<std::pair<G4IT*, double> >* FindNearest(const T* it0, const T* it);
    template<typename T> std::vector<std::pair<G4IT*, double> >* FindNearestInRange(const G4ThreeVector& pos,
            const T* it, G4double range);
    template<typename T> std::vector<std::pair<G4IT*, double> >* FindNearestInRange(const T* it0, const T* it, G4double range);
};

#endif
