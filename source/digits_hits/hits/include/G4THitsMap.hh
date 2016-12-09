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
// $Id: G4THitsMap.hh 99262 2016-09-09 13:18:19Z gcosmo $
//
#ifndef G4THitsMap_h
#define G4THitsMap_h 1

#include "G4THitsCollection.hh"
#include "globals.hh"
#include <map>

class G4StatDouble;

// class description:
//
//  This is a template class of hits map and parametrized by
// The concrete class of G4VHit. This is a uniform collection for
// a particular concrete hit class objects.
//  An intermediate layer class G4HitsMap appeared in this
// header file is used just for G4Allocator, because G4Allocator
// cannot be instansiated with a template class. Thus G4HitsMap
// class MUST NOT be directly used by the user.

template <typename T> class G4THitsMap : public G4HitsCollection 
{
  public:
      G4THitsMap();
  public: // with description
      G4THitsMap(G4String detName,G4String colNam);
      // constructor.
  public:
      virtual ~G4THitsMap();
      G4int operator==(const G4THitsMap<T> &right) const;

  public: // with description
      // Operator += between same kind of classes
      template <typename U = T,
              typename std::enable_if<std::is_same<U,T>::value,int>::type=0>
      G4THitsMap<T> & operator+=(const G4THitsMap<T> &right) const
      {
        std::map<G4int,T*> * aHitsMap = right.GetMap();
        typename std::map<G4int,T*>::iterator itr = aHitsMap->begin();
        for(; itr != aHitsMap->end(); itr++) {
	  add(itr->first, *(itr->second));
        }
        return (G4THitsMap<T>&)(*this);
      }
      // Operator for G4THitsMap<G4StatDouble> += G4THitsMap<G4double>
      template <typename U = T,
              typename std::enable_if<std::is_same<U,G4double>::value &&
                       std::is_same<T,G4StatDouble>::value,int>::type=0>
      G4THitsMap<T> & operator+=(const G4THitsMap<U> &right) const
      {
        std::map<G4int,U*> * aHitsMap = right.GetMap();
        typename std::map<G4int,U*>::iterator itr = aHitsMap->begin();
        for(; itr != aHitsMap->end(); itr++) {
          typename std::map<G4int,T*>::iterator mapItr = this->GetMap()->find(itr->first);
          if(mapItr==this->GetMap()->end())
          { (*this->GetMap())[itr->first] = new T(0.); }
	  add(itr->first, *(itr->second));
        }
        return (G4THitsMap<T>&)(*this);
      }

  public: // with description
      virtual void DrawAllHits();
      virtual void PrintAllHits();
      //  These two methods invokes Draw() and Print() methods of all of
      // hit objects stored in this map, respectively.

  public: // with description
      //  Returns a pointer to a concrete hit object.
      inline T* operator[](G4int key) const;
      //  Returns a collection map.
      inline std::map<G4int,T*>* GetMap() const
      { return (std::map<G4int,T*>*)theCollection; }

      //  Insert a hit object. Total number of hit objects stored in this
      // map is returned.
      template <typename U = T,
              typename std::enable_if<std::is_same<U,T>::value,int>::type=0>
      inline G4int add(const G4int & key, T * &aHit) const
      {
        typename std::map<G4int,T*> * theHitsMap = GetMap();
        if(theHitsMap->find(key) != theHitsMap->end()) {
          *(*theHitsMap)[key] += *aHit;
        } else {
          (*theHitsMap)[key] = aHit;
        }
        return theHitsMap->size();
      }
      template <typename U = T,
              typename std::enable_if<std::is_same<U,T>::value,int>::type=0>
      inline G4int add(const G4int & key, T &aHit) const
      {
        typename std::map<G4int,T*> * theHitsMap = GetMap();
        if(theHitsMap->find(key) == theHitsMap->end()) {
          (*theHitsMap)[key] = new T(0.);
        }
        *(*theHitsMap)[key] += aHit;
        return theHitsMap->size();
      }
      template <typename U = T,
              typename std::enable_if<std::is_same<U,G4double>::value &&
                       std::is_same<T,G4StatDouble>::value,int>::type=0>
      inline G4int add(const G4int & key, U * &aHit) const
      {
        typename std::map<G4int,T*> * theHitsMap = GetMap();
        if(theHitsMap->find(key) == theHitsMap->end()) {
          (*theHitsMap)[key] = new T(0.);
        }
        *(*theHitsMap)[key] += *aHit;
        return theHitsMap->size();
      }
      template <typename U = T,
              typename std::enable_if<std::is_same<U,G4double>::value &&
                       std::is_same<T,G4StatDouble>::value,int>::type=0>
      inline G4int add(const G4int & key, U &aHit) const
      {
        typename std::map<G4int,T*> * theHitsMap = GetMap();
        if(theHitsMap->find(key) == theHitsMap->end()) {
          (*theHitsMap)[key] = new T(0.);
        }
        *(*theHitsMap)[key] += aHit;
        return theHitsMap->size();
      }

      //  Overwrite a hit object. Total number of hit objects stored in this
      // map is returned.
      template <typename U = T,
              typename std::enable_if<std::is_same<U,T>::value,int>::type=0>
      inline G4int set(const G4int & key, T * &aHit) const
      {
        typename std::map<G4int,T*> * theHitsMap = GetMap();
        if(theHitsMap->find(key) != theHitsMap->end()) {
            delete (*theHitsMap)[key]->second;
        }
        (*theHitsMap)[key] = aHit;
        return theHitsMap->size();
      }
      template <typename U = T,
              typename std::enable_if<std::is_same<U,T>::value,int>::type=0>
      inline G4int set(const G4int & key, T &aHit) const
      {
        typename std::map<G4int,T*> * theHitsMap = GetMap();
        if(theHitsMap->find(key) == theHitsMap->end()) 
        { (*theHitsMap)[key] = new T(0.); }
        *(*theHitsMap)[key] = aHit;
        return theHitsMap->size();
      }
      template <typename U = T,
              typename std::enable_if<std::is_same<U,G4double>::value &&
                       std::is_same<T,G4StatDouble>::value,int>::type=0>
      inline G4int set(const G4int & key, U * &aHit) const
      {
        typename std::map<G4int,T*> * theHitsMap = GetMap();
        if(theHitsMap->find(key) != theHitsMap->end()) {
            delete (*theHitsMap)[key]->second;
        }
        (*theHitsMap)[key] = aHit;
        return theHitsMap->size();
      }
      template <typename U = T,
              typename std::enable_if<std::is_same<U,G4double>::value &&
                       std::is_same<T,G4StatDouble>::value,int>::type=0>
      inline G4int set(const G4int & key, U &aHit) const
      {
        typename std::map<G4int,T*> * theHitsMap = GetMap();
        if(theHitsMap->find(key) == theHitsMap->end()) 
        { (*theHitsMap)[key] = new T(0.); }
        *(*theHitsMap)[key] = aHit;
        return theHitsMap->size();
      }

      //  Returns the number of hit objects stored in this map
      inline G4int entries() const
      { return ((std::map<G4int,T*>*)theCollection)->size(); }

      // Clear the entry
      inline void clear();

  public:
    virtual G4VHit* GetHit(size_t) const {return 0;}
    virtual size_t GetSize() const
    { return ((std::map<G4int,T*>*)theCollection)->size(); }

};

template <typename T> G4THitsMap<T>::G4THitsMap()
{ 
  theCollection = (void*)new std::map<G4int,T*>;
}

template <typename T> G4THitsMap<T>::G4THitsMap(G4String detName,G4String colNam)
    : G4HitsCollection(detName,colNam)
{ 
    theCollection = (void*)new std::map<G4int,T*>;
}

template <typename T> G4THitsMap<T>::~G4THitsMap()
{
  typename std::map<G4int,T*> * theHitsMap = GetMap();
  typename std::map<G4int,T*>::iterator itr = theHitsMap->begin();
  for(; itr != theHitsMap->end(); itr++) {
      delete itr->second;
  }

  delete theHitsMap;
}

template <typename T> G4int G4THitsMap<T>::operator==(const G4THitsMap<T> &right) const
{ return (collectionName==right.collectionName); }

template <typename T> inline T* 
G4THitsMap<T>::operator[](G4int key) const {
    std::map<G4int,T*> * theHitsMap = GetMap();
    if(theHitsMap->find(key) != theHitsMap->end()) {
	return theHitsMap->find(key)->second;
    } else {
	return 0;
    }
}

template <typename T> void G4THitsMap<T>::DrawAllHits() 
{;}

template <typename T> void G4THitsMap<T>::PrintAllHits() 
{
 G4cout << "G4THitsMap " << SDname << " / " << collectionName << " --- " << entries() << " entries" << G4endl;
/*----- commented out for the use-case where <T> cannot be initialized
        to be zero or does not support += operator.
 std::map<G4int,T*> * theHitsMap = GetMap();
 typename std::map<G4int, T*>::iterator itr = theHitsMap->begin();
 T sum = 0.;
 for(; itr != theHitsMap->end(); itr++) {
  ///////////////////////////////G4cout << "  " << itr->first << " : " << *(itr->second) << G4endl;
  sum += *(itr->second);
 }
 G4cout << "             Total : " << sum << G4endl;
----------------------------------------------------------------------*/
}

template <typename T> void G4THitsMap<T>::clear() {

    std::map<G4int,T*> * theHitsMap = GetMap();
    typename std::map<G4int, T*>::iterator itr = theHitsMap->begin();
    for(; itr != theHitsMap->end(); itr++) {
	delete itr->second;
    }
    theHitsMap->clear();

}

#endif

