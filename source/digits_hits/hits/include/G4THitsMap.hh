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
//
//
// $Id: G4THitsMap.hh,v 1.5 2005/11/18 22:47:06 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
#ifndef G4THitsMap_h
#define G4THitsMap_h 1

#include "G4THitsCollection.hh"
#include "globals.hh"
#include <map>

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
      G4THitsMap<T> & operator+=(const G4THitsMap<T> &right) const;

  public: // with description
      virtual void DrawAllHits();
      virtual void PrintAllHits();
      //  These two methods invokes Draw() and Print() methods of all of
      // hit objects stored in this map, respectively.

  public: // with description
      inline T* operator[](G4int key) const;

      //  Returns a pointer to a concrete hit object.
      inline std::map<G4int,T*>* GetMap() const
      { return (std::map<G4int,T*>*)theCollection; }
      //  Returns a collection map.
      inline G4int add(const G4int & key, T * &aHit) const;
      inline G4int add(const G4int & key, T &aHit) const;
      //  Insert a hit object. Total number of hit objects stored in this
      // map is returned.
      inline G4int set(const G4int & key, T * &aHit) const;
      inline G4int set(const G4int & key, T &aHit) const;
      //  Overwrite a hit object. Total number of hit objects stored in this
      // map is returned.
      inline G4int entries() const
      { return ((std::map<G4int,T*>*)theCollection)->size(); }
      //  Returns the number of hit objects stored in this map
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

template <typename T> G4THitsMap<T> &
G4THitsMap<T>::operator+=(const G4THitsMap<T> &right) const
{
    std::map<G4int,T*> * aHitsMap = right.GetMap();
    typename std::map<G4int,T*>::iterator itr = aHitsMap->begin();
    for(; itr != aHitsMap->end(); itr++) {
	add(itr->first, *(itr->second));
    }
    return (G4THitsMap<T>&)(*this);
}

template <typename T> inline T* 
G4THitsMap<T>::operator[](G4int key) const {
    std::map<G4int,T*> * theHitsMap = GetMap();
    if(theHitsMap->find(key) != theHitsMap->end()) {
	return theHitsMap->find(key)->second;
    } else {
	return 0;
    }
}

template <typename T> inline G4int
G4THitsMap<T>::add(const G4int & key, T * &aHit) const {

    typename std::map<G4int,T*> * theHitsMap = GetMap();
    if(theHitsMap->find(key) != theHitsMap->end()) {
	*(*theHitsMap)[key] += *aHit;
    } else {
	(*theHitsMap)[key] = aHit;
    }
    return theHitsMap->size();
}

template <typename T> inline G4int
G4THitsMap<T>::add(const G4int & key, T &aHit) const {

    typename std::map<G4int,T*> * theHitsMap = GetMap();
    if(theHitsMap->find(key) != theHitsMap->end()) {
	*(*theHitsMap)[key] += aHit;
    } else {
	T * hit = new T;
	*hit = aHit;
	(*theHitsMap)[key] = hit;
    }

    return theHitsMap->size();
}

template <typename T> inline G4int
G4THitsMap<T>::set(const G4int & key, T * &aHit) const {
                                                                                             
    typename std::map<G4int,T*> * theHitsMap = GetMap();
    if(theHitsMap->find(key) != theHitsMap->end()) {
        delete (*theHitsMap)[key]->second;
    }
    (*theHitsMap)[key] = aHit;
    return theHitsMap->size();
}
                                                                                             
template <typename T> inline G4int
G4THitsMap<T>::set(const G4int & key, T &aHit) const {
                                                                                             
    typename std::map<G4int,T*> * theHitsMap = GetMap();
    if(theHitsMap->find(key) != theHitsMap->end()) {
        *(*theHitsMap)[key] = aHit;
    } else {
        T * hit = new T;
        *hit = aHit;
        (*theHitsMap)[key] = hit;
    }
                                                                                             
    return theHitsMap->size();
}
                                                                                             
template <typename T> void G4THitsMap<T>::DrawAllHits() 
{;}

template <typename T> void G4THitsMap<T>::PrintAllHits() 
{;}

template <typename T> void G4THitsMap<T>::clear() {

    std::map<G4int,T*> * theHitsMap = GetMap();
    typename std::map<G4int, T*>::iterator itr = theHitsMap->begin();
    for(; itr != theHitsMap->end(); itr++) {
	delete itr->second;
    }
    theHitsMap->clear();

}

#endif

