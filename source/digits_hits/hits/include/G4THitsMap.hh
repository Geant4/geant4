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
// $Id: G4THitsMap.hh 103976 2017-05-05 12:22:53Z gcosmo $
//
#ifndef G4THitsMap_h
#define G4THitsMap_h 1

#include "G4THitsCollection.hh"
#include "globals.hh"

#include <map>
#include <unordered_map>

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

template <typename T, typename Map_t = std::map<G4int, T*> >
class G4VTHitsMap : public G4HitsCollection
{
private:
    typedef std::multimap<G4int, T*>            mmap_t;
    typedef std::pair<G4int, T*>                pair_t;
    typedef std::unordered_multimap<G4int, T*>  uommap_t;

    #define is_same_t(_Tp, _Up) std::is_same<_Tp, _Up>::value
    #define is_multimap_t(_Mp) std::is_same<_Mp, mmap_t>::value
    #define is_uommap_t(_Mp) std::is_same<_Mp, uommap_t>::value
    #define is_mmap_t(_Mp) (is_multimap_t(_Mp) || is_uommap_t(_Mp))
    #define is_fundamental_t(_Tp) std::is_fundamental<_Tp>::value

    template <bool _Bp, typename _Tp = void>
    using enable_if_t = typename std::enable_if<_Bp, _Tp>::type;

    // ensure fundamental types are initialized to zero
    template <typename U = T, enable_if_t<is_fundamental_t(U), int> = 0>
    T* allocate() const
    { return new T(0.); }
    // non-fundamental types should set values to appropriate values
    // and avoid issues such as:
    //   G4StatDouble stat(0.); stat += 1.0; gives n == 2;
    template <typename U = T, enable_if_t<! is_fundamental_t(U), int> = 0>
    T* allocate() const
    { return new T(); }

public:
    typedef G4VTHitsMap<T, Map_t>               this_type;
    typedef T                                   value_type;
    typedef Map_t                               map_type;
    typedef typename map_type::iterator         iterator;
    typedef typename map_type::const_iterator   const_iterator;

public: // with description
    // generic constructor
    G4VTHitsMap();
    // det + collection description constructor
    G4VTHitsMap(G4String detName, G4String colNam);
    // destructor
    virtual ~G4VTHitsMap();
    // equivalence operator
    G4int operator==(const G4VTHitsMap<T, Map_t> &right) const;

    //------------------------------------------------------------------------//
    // Generic operator += where add(...) overloads handle various
    //  U and MapU_t types
    //------------------------------------------------------------------------//
    template <typename U, typename MapU_t>
    this_type& operator+=(const G4VTHitsMap<U, MapU_t>& right) const
    {
        MapU_t* aHitsMap = right.GetMap();
        for(auto itr = aHitsMap->begin(); itr != aHitsMap->end(); itr++)
            add<U, map_type>(itr->first, *(itr->second));
        return (this_type&)(*this);
    }
    //------------------------------------------------------------------------//

public: // with description
    virtual void DrawAllHits();
    virtual void PrintAllHits();
    //  These two methods invokes Draw() and Print() methods of all of
    //  hit objects stored in this map, respectively.

public: // with description
    //  Returns a pointer to a concrete hit object.
    inline Map_t* GetMap() const
    { return (Map_t*)theCollection; }
    //  Overwrite a hit object. Total number of hit objects stored in this
    // map is returned.
    inline G4int entries() const
    { return ((Map_t*)theCollection)->size(); }
    //  Returns the number of hit objects stored in this map
    inline void clear();

    iterator begin()              { return GetMap()->begin(); }
    iterator end()                { return GetMap()->end(); }
    const_iterator begin() const  { return GetMap()->begin(); }
    const_iterator end() const    { return GetMap()->end(); }
    const_iterator cbegin() const { return GetMap()->cbegin(); }
    const_iterator cend() const   { return GetMap()->cend(); }

public:
    virtual G4VHit* GetHit(size_t) const {return 0;}
    virtual size_t GetSize() const
    { return ((Map_t*)theCollection)->size(); }

public:
    //------------------------------------------------------------------------//
    //  Add/Insert a hit object. Total number of hit objects stored in this
    //  map is returned.
    //------------------------------------------------------------------------//
    //  Standard map overload for same type
    //------------------------------------------------------------------------//
    // here we don't use allocate() since instances like G4Colour() == white
    // and += adds to white (not correct)
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<is_same_t(U, T) && ! is_mmap_t(MapU_t), int> = 0>
    G4int add(const G4int& key, U*& aHit) const
    {
        map_type* theHitsMap = GetMap();
        if(theHitsMap->find(key) == theHitsMap->end())
            theHitsMap->insert(pair_t(key, new T(*aHit)));
        else
            *theHitsMap->find(key)->second += *aHit;
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Multimap overload for same type T
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(is_same_t(U, T) && is_mmap_t(MapU_t)), int> = 0>
    G4int add(const G4int& key, U*& aHit) const
    {
        map_type* theHitsMap = GetMap();
        theHitsMap->insert(pair_t(key, aHit));
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Multimap overload for different types
    //      assumes type T has overload of += operator for U
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(!is_same_t(U, T) && is_mmap_t(MapU_t)), int> = 0>
    G4int add(const G4int& key, U*& aHit) const
    {
        map_type* theHitsMap = GetMap();
        T* hit = allocate();
        *hit += *aHit;
        theHitsMap->insert(pair_t(key, hit));
        return theHitsMap->size();
    }

public:
    //------------------------------------------------------------------------//
    //  Standard map overload for same type
    //      assumes type T has overload of += operator for U
    //------------------------------------------------------------------------//
    // here we don't use allocate() since instances like G4Colour() == white
    // and += adds to white (not correct)
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(is_same_t(U, T) && ! is_mmap_t(MapU_t)), int> = 0>
    G4int add(const G4int& key, U& aHit) const
    {
        map_type* theHitsMap = GetMap();
        if(theHitsMap->find(key) == theHitsMap->end())
            theHitsMap->insert(pair_t(key, new T(aHit)));
        else
            *theHitsMap->find(key)->second += aHit;
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Standard map overload for different type
    //      assumes type T has overload of += operator for U
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(! is_same_t(U, T) && ! is_mmap_t(MapU_t)), int> = 0>
    G4int add(const G4int& key, U& aHit) const
    {
        map_type* theHitsMap = GetMap();
        if(theHitsMap->find(key) == theHitsMap->end())
            theHitsMap->insert(pair_t(key, allocate()));
        *theHitsMap->find(key)->second += aHit;
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Multimap overload for same type T
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(is_same_t(U, T) && is_mmap_t(MapU_t)), int> = 0>
    G4int add(const G4int& key, U& aHit) const
    {
        map_type* theHitsMap = GetMap();
        theHitsMap->insert(pair_t(key, new T(aHit)));
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Multimap overload for different types
    //      assumes type T has overload of += operator for U
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(!is_same_t(U, T) && is_mmap_t(MapU_t)), int> = 0>
    G4int add(const G4int& key, U& aHit) const
    {
        map_type* theHitsMap = GetMap();
        T* hit = allocate();
        *hit += aHit;
        theHitsMap->insert(pair_t(key, hit));
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//

public:
    //------------------------------------------------------------------------//
    //  Set a hit object. Total number of hit objects stored in this
    //  map is returned.
    //------------------------------------------------------------------------//
    //  Standard overload for same type T
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(is_same_t(U, T) && ! is_mmap_t(MapU_t)), int> = 0>
    inline G4int set(const G4int& key, U*& aHit) const
    {
        map_type* theHitsMap = GetMap();
        if(theHitsMap->find(key) != theHitsMap->end())
            delete theHitsMap->find(key)->second;
        theHitsMap->find(key)->second = aHit;
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Multimap overload for same type T
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(is_same_t(U, T) && is_mmap_t(MapU_t)), int> = 0>
    inline G4int set(const G4int& key, U*& aHit) const
    {
        map_type* theHitsMap = GetMap();
        if(theHitsMap->find(key) != theHitsMap->end())
            theHitsMap->insert(pair_t(key, aHit));
        else
        {
            delete theHitsMap->find(key)->second;
            theHitsMap->find(key)->second = aHit;
        }
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Standard map overload for different types
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(!is_same_t(U, T) && ! is_mmap_t(MapU_t)), int> = 0>
    inline G4int set(const G4int & key, U*& aHit) const
    {
        map_type* theHitsMap = GetMap();
        T* hit = nullptr;
        if(theHitsMap->find(key) == theHitsMap->end())
            theHitsMap->insert(std::make_pair(key, hit = allocate()));
        else
            hit = theHitsMap->find(key)->second;
        *hit += *aHit;
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Multimap overload for different types
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(!is_same_t(U, T) && is_mmap_t(MapU_t)), int> = 0>
    inline G4int set(const G4int& key, U*& aHit) const
    {
        map_type* theHitsMap = GetMap();
        T* hit = allocate();
        *hit += *aHit;
        if(theHitsMap->find(key) != theHitsMap->end())
            theHitsMap->insert(pair_t(key, hit));
        else
        {
            delete theHitsMap->find(key)->second;
            theHitsMap->find(key)->second = hit;
        }
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//

public:
    //------------------------------------------------------------------------//
    //  Set a hit object. Total number of hit objects stored in this
    //  map is returned.
    //------------------------------------------------------------------------//
    //  Standard overload for same type T
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(is_same_t(U, T) && ! is_mmap_t(MapU_t)), int> = 0>
    inline G4int set(const G4int& key, U& aHit) const
    {
        map_type* theHitsMap = GetMap();
        T* hit = nullptr;
        if(theHitsMap->find(key) != theHitsMap->end())
            hit = theHitsMap->find(key)->second;
        else
            theHitsMap->insert(pair_t(key, hit = allocate()));
        *hit = aHit;
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Multimap overload for same type T
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(is_same_t(U, T) && is_mmap_t(MapU_t)), int> = 0>
    inline G4int set(const G4int& key, U& aHit) const
    {
        map_type* theHitsMap = GetMap();
        if(theHitsMap->find(key) != theHitsMap->end())
            *theHitsMap->find(key)->second = aHit;
        else
            theHitsMap->insert(pair_t(key, new T(aHit)));
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Standard map overload for different types
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(!is_same_t(U, T) && ! is_mmap_t(MapU_t)), int> = 0>
    inline G4int set(const G4int & key, U& aHit) const
    {
        map_type* theHitsMap = GetMap();
        T* hit = nullptr;
        if(theHitsMap->find(key) == theHitsMap->end())
            theHitsMap->insert(std::make_pair(key, hit = allocate()));
        else
            hit = theHitsMap->find(key)->second;
        *hit += aHit;
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//
    //  Multimap overload for different types
    //------------------------------------------------------------------------//
    template <typename U = T,
              typename MapU_t = Map_t,
              enable_if_t<(!is_same_t(U, T) && is_mmap_t(MapU_t)), int> = 0>
    inline G4int set(const G4int& key, U& aHit) const
    {
        map_type* theHitsMap = GetMap();
        T* hit = allocate();
        *hit += aHit;
        if(theHitsMap->find(key) != theHitsMap->end())
            *theHitsMap->find(key)->second = *hit;
        else
            theHitsMap->insert(pair_t(key, hit));
        return theHitsMap->size();
    }
    //------------------------------------------------------------------------//

public:
    //------------------------------------------------------------------------//
    //  Enable bracket operator. Return pointer to data indexed by key or
    //      last occurring instance of pointer to data index by key in the
    //      case of a multimap
    //------------------------------------------------------------------------//
    template <typename MapU_t = Map_t,
              enable_if_t<! is_mmap_t(MapU_t), int> = 0>
    T* operator[](G4int key) const
    {
        map_type* theHitsMap = GetMap();
        if(theHitsMap->find(key) != theHitsMap->end())
            return theHitsMap->find(key)->second;
        return nullptr;
    }
    //------------------------------------------------------------------------//
    template <typename MapU_t = Map_t,
              enable_if_t<is_mmap_t(MapU_t), int> = 0>
    T* operator[](G4int key) const
    {
#ifdef G4VERBOSE
        static bool _first = true;
        if(_first)
        {
            _first = false;
            G4Exception("G4THitsMap operator[]", "calling [] on multimap",
                        JustWarning, "Returning the last matching entry");
        }
#endif
        map_type* theHitsMap = GetMap();
        iterator itr = theHitsMap->find(key);
        if(itr != theHitsMap->end())
        {
            std::advance(itr, theHitsMap->count(key)-1);
            return itr->second;
        }
        return nullptr;
    }
    //------------------------------------------------------------------------//

    #undef is_same_t
    #undef is_multimap_t
    #undef is_uommap_t
    #undef is_mmap_t
    #undef is_fundamental_t
};

//============================================================================//

template <typename T, typename Map_t>
G4VTHitsMap<T, Map_t>::G4VTHitsMap()
{
  theCollection = (void*)new Map_t;
}

//============================================================================//

template <typename T, typename Map_t>
G4VTHitsMap<T, Map_t>::G4VTHitsMap(G4String detName,G4String colNam)
    : G4HitsCollection(detName,colNam)
{
    theCollection = (void*)new Map_t;
}

//============================================================================//

template <typename T, typename Map_t>
G4VTHitsMap<T, Map_t>::~G4VTHitsMap()
{
  map_type* theHitsMap = GetMap();
  for(iterator itr = theHitsMap->begin(); itr != theHitsMap->end(); itr++)
      delete itr->second;
  delete theHitsMap;
}

//============================================================================//

template <typename T, typename Map_t>
G4int G4VTHitsMap<T, Map_t>::operator==(const G4VTHitsMap<T, Map_t> &right) const
{
    return (collectionName==right.collectionName);
}

//============================================================================//

template <typename T, typename Map_t>
void G4VTHitsMap<T, Map_t>::DrawAllHits()
{;}

//============================================================================//

template <typename T, typename Map_t>
void G4VTHitsMap<T, Map_t>::PrintAllHits()
{
    G4cout << "G4THitsMap " << SDname << " / " << collectionName << " --- "
           << entries() << " entries" << G4endl;
    /*----- commented out for the use-case where <T> cannot be initialized
            to be zero or does not support += operator.
     Map_t * theHitsMap = GetMap();
     typename Map_t::iterator itr = theHitsMap->begin();
     T sum = 0.;
     for(; itr != theHitsMap->end(); itr++) {
      G4cout << "  " << itr->first << " : "
             << *(itr->second) << G4endl;
      sum += *(itr->second);
     }
     G4cout << "             Total : " << sum << G4endl;
    ----------------------------------------------------------------------*/
}

//============================================================================//

template <typename T, typename Map_t>
void G4VTHitsMap<T, Map_t>::clear()
{
    Map_t * theHitsMap = GetMap();
    for(iterator itr = theHitsMap->begin(); itr != theHitsMap->end(); itr++)
        delete itr->second;
    theHitsMap->clear();
}

//============================================================================//
//                                                                            //
//                                                                            //
//                      Helpers for different map types                       //
//                                                                            //
//                                                                            //
//============================================================================//

template <typename _Tp>
class G4THitsMap : public G4VTHitsMap<_Tp, std::map<G4int, _Tp*>>
{
public:
    typedef G4VTHitsMap<_Tp, std::map<G4int, _Tp*>> parent_type;

public:
    G4THitsMap() : parent_type() { }
    G4THitsMap(G4String detName,G4String colName)
    : parent_type(detName, colName) { }

    using parent_type::operator +=;
    using parent_type::operator ==;
    using parent_type::operator [];
    using parent_type::DrawAllHits;
    using parent_type::PrintAllHits;
    using parent_type::GetMap;
    using parent_type::entries;
    using parent_type::clear;
    using parent_type::begin;
    using parent_type::end;
    using parent_type::cbegin;
    using parent_type::cend;
    using parent_type::GetHit;
    using parent_type::GetSize;
    using parent_type::add;
    using parent_type::set;
};

//============================================================================//

template <typename _Tp>
class G4THitsMultiMap : public G4VTHitsMap<_Tp, std::multimap<G4int, _Tp*>>
{
public:
    typedef G4VTHitsMap<_Tp, std::multimap<G4int, _Tp*>> parent_type;

public:
    G4THitsMultiMap() : parent_type() { }
    G4THitsMultiMap(G4String detName, G4String colName)
    : parent_type(detName, colName) { }

    using parent_type::operator +=;
    using parent_type::operator ==;
    using parent_type::operator [];
    using parent_type::DrawAllHits;
    using parent_type::PrintAllHits;
    using parent_type::GetMap;
    using parent_type::entries;
    using parent_type::clear;
    using parent_type::begin;
    using parent_type::end;
    using parent_type::cbegin;
    using parent_type::cend;
    using parent_type::GetHit;
    using parent_type::GetSize;
    using parent_type::add;
    using parent_type::set;
};

//============================================================================//

template <typename _Tp>
class G4THitsUnorderedMap
        : public G4VTHitsMap<_Tp, std::unordered_map<G4int, _Tp*>>
{
public:
    typedef G4VTHitsMap<_Tp, std::unordered_map<G4int, _Tp*>> parent_type;

public:
    G4THitsUnorderedMap() : parent_type() { }
    G4THitsUnorderedMap(G4String detName, G4String colName)
    : parent_type(detName, colName) { }

    using parent_type::operator +=;
    using parent_type::operator ==;
    using parent_type::operator [];
    using parent_type::DrawAllHits;
    using parent_type::PrintAllHits;
    using parent_type::GetMap;
    using parent_type::entries;
    using parent_type::clear;
    using parent_type::begin;
    using parent_type::end;
    using parent_type::cbegin;
    using parent_type::cend;
    using parent_type::GetHit;
    using parent_type::GetSize;
    using parent_type::add;
    using parent_type::set;
};

//============================================================================//

template <typename _Tp>
class G4THitsUnorderedMultiMap
        : public G4VTHitsMap<_Tp, std::unordered_multimap<G4int, _Tp*>>
{
public:
    typedef G4VTHitsMap<_Tp, std::unordered_multimap<G4int, _Tp*>> parent_type;

public:
    G4THitsUnorderedMultiMap() : parent_type() { }
    G4THitsUnorderedMultiMap(G4String detName,G4String colName)
    : parent_type(detName, colName) { }

    using parent_type::operator +=;
    using parent_type::operator ==;
    using parent_type::operator [];
    using parent_type::DrawAllHits;
    using parent_type::PrintAllHits;
    using parent_type::GetMap;
    using parent_type::entries;
    using parent_type::clear;
    using parent_type::begin;
    using parent_type::end;
    using parent_type::cbegin;
    using parent_type::cend;
    using parent_type::GetHit;
    using parent_type::GetSize;
    using parent_type::add;
    using parent_type::set;
};

//============================================================================//

#endif
