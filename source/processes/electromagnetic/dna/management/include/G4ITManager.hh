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
// $Id: G4ITManager.hh 65022 2012-11-12 16:43:12Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4ITManager_hh
#define G4ITManager_hh 1

#include "globals.hh"
#include <map>
#include "G4AllITManager.hh"
#include "G4ITBox.hh"
#include "G4KDTree.hh"
#include "G4Track.hh"

/**
  * G4VITManager is just a virtual interface for G4ITManager.
  * For more details, please have a look at the description
  * of ITManager.
  */

class G4VITManager
{    
protected :
    G4ITType fType ;
    G4int fVerbose;

public :
    G4VITManager();
    virtual ~G4VITManager(){;}

    void SetVerboseLevel(G4int level)
    {
        fVerbose = level;
    }
    G4int GetVerboseLevel()
    {
        return fVerbose;
    }

    virtual void UpdatePositionMap() = 0;
    virtual void CreateTree() {;}

    virtual void Push(G4Track*) = 0 ;
    G4ITType GetITType()
    {
        return fType;
    }

    G4ITBox* GetBox(const G4Track* track)
    {
        return GetBox(GetIT(track));
    }

    virtual G4ITBox* GetBox(const G4IT*) = 0;

    // Navigate between boxes
    virtual G4ITBox* GetFirstBox() = 0;
    virtual G4ITBox* GetNextBox(G4ITBox*) = 0;
    virtual G4ITBox* GetLastBox() = 0;

public :

    class iterator
    {
    public :
        iterator(G4ITBox*);
        virtual ~iterator(){;}
        virtual G4bool begin();
        virtual G4bool end();
        iterator& operator= (const iterator& i);
        iterator& operator++(G4int);
        G4IT* operator*();
        G4ITBox* GetBox();

    protected :
        //_____________________________________
        // Print "G4IT* friendIT" status :
        void PrintNext() const;
        //_____________________________________
        // Attributes
        G4IT* fNextIT; // the one you are looking reactants for
        G4ITBox* fCurrentBox;
    };

    class allbox_iterator : public iterator
    {
    public :
        allbox_iterator(G4ITType);
        allbox_iterator(G4VITManager*);
        virtual ~allbox_iterator() {;}
        allbox_iterator& operator= (const allbox_iterator& i);
        allbox_iterator& operator++(G4int);

    protected :
        G4VITManager* fInstance ;

    };

    class const_iterator : public G4VITManager::iterator
    {
    public :
        const_iterator(G4ITBox*);
        virtual ~const_iterator(){;}
        const G4IT* operator*();
    };
};

/**
  * G4ITManager is able to save into different boxes
  * the ITs that will be used in the simulation.
  * It is a stack-like.
  * Those boxes are used to fill a tree which helps
  * finding the closest neighboor.
  */

template<typename T>
class G4ITManager : public G4VITManager
{
    static G4ITManager<T> *      fInstance;
    G4ITManager<T>();

    typedef std::map<T,G4ITBox* > BoxMap;
    BoxMap fBox;

    typedef std::map<T, G4KDTree* > TreeMap;
    TreeMap fTree;
    TreeMap fPrevious_tree;

public :
    static G4ITManager<T> *      Instance();
    virtual ~G4ITManager();
    virtual void Push(G4Track*);

    //Make a method to extract number of IT from Box
    G4int NbElements(const G4IT*);

    void EraseABox(T*);
    void EraseABox(G4ITBox*);

    void SetVerboseLevel(G4int level)
    {
        fVerbose = level;
    }
    G4int GetVerboseLevel()
    {
        return fVerbose;
    }

    virtual void UpdatePositionMap();
    static void iUpdatePositionMap();

    G4KDTreeResultHandle FindNearestInRange(const T*, const T*, G4double);
    G4KDTreeResultHandle FindNearest(const G4ThreeVector&, const T* it);
    G4KDTreeResultHandle FindNearest(const T* it0, const T* it);
    G4KDTreeResultHandle FindNearestInRange(const G4ThreeVector&, const T*, G4double);

    inline G4ITBox* GetBox(const T* IT)
    {
        typename BoxMap::const_iterator it = fBox.find(*IT);
        if(it == fBox.end()) return 0;
        return it->second;
    }

    inline virtual G4ITBox* GetBox(const G4IT* IT)
    {
        const T* myIT = dynamic_cast<const T*>(IT);

        if(myIT == 0)
        {
            G4ExceptionDescription exceptionDescription ("You are requested a bad IT");
            G4Exception("G4ITManager::GetBox","ITManager001",
                        FatalErrorInArgument,exceptionDescription);
            return 0; // coverity
        }

        return GetBox(myIT);
    }

    virtual G4ITBox* GetFirstBox()
    {
        typename BoxMap::iterator it = fBox.begin();
        if(it != fBox.end())
        {
            return it->second;
        }
        return 0;
    }

    virtual G4ITBox* GetNextBox(G4ITBox* box)
    {
        if(box)
        {
            return box->GetNextBox();
        }
        return 0;
    }

    virtual G4ITBox* GetLastBox()
    {
        typename BoxMap::reverse_iterator it = fBox.rbegin();
        if(it != fBox.rend())
        {
            return it->second;
        }
        return 0;
    }

};

#ifdef TEMPLATE
#undef TEMPLATE
#endif

#define TEMPLATE template<typename T>
#define G4ITMANAGER G4ITManager<T>

#include "G4ITManager.icc"

#undef TEMPLATE
#undef G4ITMANAGER

#endif
