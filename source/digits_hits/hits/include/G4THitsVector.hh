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
//
#ifndef G4THitsVector_h
#define G4THitsVector_h 1

#include "G4THitsCollection.hh"
#include "G4THitsMap.hh"
#include "globals.hh"

#include <deque>
#include <map>
#include <unordered_map>
#include <vector>

// class description:
//
//  This is a template class of hits Vector and parametrized by
// The concrete class of G4VHit. This is a uniform collection for
// a particular concrete hit class objects.
//  An intermediate layer class G4HitsVector appeared in this
// header file is used just for G4Allocator, because G4Allocator
// cannot be instansiated with a template class. Thus G4HitsVector
// class MUST NOT be directly used by the user.

template <typename T, typename Vector_t = std::deque<T*>>
class G4VTHitsVector : public G4HitsCollection
{
 public:
  using this_type = G4VTHitsVector<T, Vector_t>;
  using value_type = T;
  using vector_type = Vector_t;
  using iterator = typename vector_type::iterator;
  using const_iterator = typename vector_type::const_iterator;

  using store_type = typename Vector_t::value_type;
  using pair_t = std::pair<G4int, store_type>;
  using map_t = std::map<G4int, store_type>;
  using uomap_t = std::unordered_map<G4int, store_type>;
  using mmap_t = std::multimap<G4int, store_type>;
  using uommap_t = std::unordered_multimap<G4int, store_type>;

 private:
#define is_same_t(_Tp, _Up) std::is_same<_Tp, _Up>::value
#define is_fundamental_t(_Tp) std::is_fundamental<_Tp>::value

#define is_std_map_t(_Mp) std::is_same<_Mp, map_t>::value
#define is_std_uomap_t(_Mp) std::is_same<_Mp, uomap_t>::value
#define is_std_mmap_t(_Mp) std::is_same<_Mp, mmap_t>::value
#define is_std_uommap_t(_Mp) std::is_same<_Mp, uommap_t>::value
#define is_map_t(_Mp) \
  (is_std_map_t(_Mp) ||\ is_std_mmap_t(_Mp) || \ is_std_uomap_t(_Mp) || \ is_std_uommap_t(_Mp))
#define is_pointer_t(_Tp) std::is_pointer<_Tp>::value
#define scast(_Tp) static_cast<_Tp>

  template <bool _Bp, typename _Tp = void>
  using enable_if_t = typename std::enable_if<_Bp, _Tp>::type;

 public:
  // generic constructor
  G4VTHitsVector(G4int init_size = 0);
  // det + collection description constructor
  G4VTHitsVector(G4String detName, G4String colNam, G4int init_size = 0);
  // destructor
  ~G4VTHitsVector() override;
  // equivalence operator
  G4bool operator==(const this_type& rhs) const;

  void DrawAllHits() override;
  void PrintAllHits() override;
  //  These two methods invokes Draw() and Print() methods of all of
  // hit objects stored in this map, respectively.

  // Generic iteration
  inline Vector_t* GetContainer() const { return scast(Vector_t*)(theCollection); }

  inline typename Vector_t::size_type size() { return GetContainer()->size(); }

  inline typename Vector_t::size_type GetIndex(iterator itr) { return std::distance(begin(), itr); }

  inline typename Vector_t::size_type GetIndex(const_iterator itr) const
  {
    return std::distance(begin(), itr);
  }

  template <typename U = store_type, enable_if_t<(is_pointer_t(U)), G4int> = 0>
  inline T* GetObject(G4int idx) const
  {
    return (idx < (G4int)GetContainer()->size()) ? (*GetContainer())[idx] : nullptr;
  }

  template <typename U = store_type, enable_if_t<(is_pointer_t(U)), G4int> = 0>
  inline T* GetObject(iterator itr) const
  {
    return (*itr);
  }

  template <typename U = store_type, enable_if_t<(is_pointer_t(U)), G4int> = 0>
  inline const T* GetObject(const_iterator itr) const
  {
    return (*itr);
  }

  template <typename U = store_type, enable_if_t<(! is_pointer_t(U)), G4int> = 0>
  inline T* GetObject(G4int idx) const
  {
    return (idx < (G4int)GetContainer()->size()) ? &(*GetContainer())[idx] : nullptr;
  }

  template <typename U = store_type, enable_if_t<(! is_pointer_t(U)), G4int> = 0>
  inline T* GetObject(iterator itr) const
  {
    return &(*itr);
  }

  template <typename U = store_type, enable_if_t<(! is_pointer_t(U)), G4int> = 0>
  inline const T* GetObject(const_iterator itr) const
  {
    return &(*itr);
  }

  iterator begin() { return GetContainer()->begin(); }
  iterator end() { return GetContainer()->end(); }
  const_iterator begin() const { return GetContainer()->begin(); }
  const_iterator end() const { return GetContainer()->end(); }
  const_iterator cbegin() const { return GetContainer()->cbegin(); }
  const_iterator cend() const { return GetContainer()->cend(); }

  //  Returns a pointer to a concrete hit object.
  inline Vector_t* GetVector() const { return scast(Vector_t*)(theCollection); }

  //  Overwrite a hit object. Total number of hit objects stored in this
  // map is returned.
  inline std::size_t entries() const { return (scast(Vector_t*)(theCollection))->size(); }

  //  Returns the number of hit objects stored in this map
  inline void clear();

  G4VHit* GetHit(std::size_t) const override { return nullptr; }
  std::size_t GetSize() const override { return (scast(Vector_t*)(theCollection))->size(); }

  inline map_t* GetMap() const;

 public:
  //------------------------------------------------------------------------//
  //  POINTER TYPE
  //------------------------------------------------------------------------//
  // ensure fundamental types are initialized to zero
  template <typename U = T, typename V = store_type,
    enable_if_t<(is_fundamental_t(U) && is_pointer_t(V)), G4int> = 0>
  store_type allocate() const
  {
    return new T(0.);
  }

  // non-fundamental types should set values to appropriate values
  // and G4StatDouble stat(0.); stat += 1.0; gives n == 2;
  template <typename U = T, typename V = store_type,
    enable_if_t<(! is_fundamental_t(U) && is_pointer_t(V)), G4int> = 0>
  store_type allocate() const
  {
    return new T();
  }

  // ensure fundamental types are initialized to zero
  template <typename U = store_type, enable_if_t<(is_pointer_t(U)), G4int> = 0>
  store_type null() const
  {
    return nullptr;
  }

  //------------------------------------------------------------------------//
  //  NON-POINTER TYPE
  //------------------------------------------------------------------------//
  // ensure fundamental types are initialized to zero
  template <typename U = T, typename V = store_type,
    enable_if_t<(is_fundamental_t(U) && ! is_pointer_t(V)), G4int> = 0>
  store_type allocate() const
  {
    return T(0.);
  }
  // non-fundamental types should set values to appropriate values
  // and G4StatDouble stat(0.); stat += 1.0; gives n == 2;
  template <typename U = T, typename V = store_type,
    enable_if_t<(! is_fundamental_t(U) && ! is_pointer_t(V)), G4int> = 0>
  store_type allocate() const
  {
    return T();
  }

  // ensure fundamental types are initialized to zero
  template <typename U = store_type, enable_if_t<(! is_pointer_t(U)), G4int> = 0>
  store_type null() const
  {
    return store_type();
  }

 public:
  //------------------------------------------------------------------------//
  // Generic operator += where add(...) overloads handle various
  //  U and VectorU_t types
  //------------------------------------------------------------------------//
  template <typename U, typename VectorU_t,
    enable_if_t<(is_pointer_t(typename VectorU_t::value_type)), G4int> = 0>
  this_type& operator+=(const G4VTHitsVector<U, VectorU_t>& right) const
  {
    VectorU_t* aHitsVector = right.GetVector();
    for (auto itr = aHitsVector->begin(); itr != aHitsVector->end(); ++itr) {
      auto _ptr = (*itr) ? (*itr) : null();
      if (_ptr) add<U>(std::distance(aHitsVector->begin(), itr), *_ptr);
    }
    return static_cast<this_type&>(*(const_cast<this_type*>(this)));
  }
  //------------------------------------------------------------------------//
  template <typename U, typename VectorU_t,
    enable_if_t<(! is_pointer_t(typename VectorU_t::value_type)), G4int> = 0>
  this_type& operator+=(const G4VTHitsVector<U, VectorU_t>& right) const
  {
    VectorU_t* aHitsVector = right.GetVector();
    for (auto itr = aHitsVector->begin(); itr != aHitsVector->end(); ++itr) {
      auto _ptr = (*itr) ? (*itr) : allocate();
      add<U>(std::distance(aHitsVector->begin(), itr), _ptr);
    }
    return static_cast<this_type&>(*(const_cast<this_type*>(this)));
  }
  //------------------------------------------------------------------------//
  template <typename U, typename MapU_t,
    enable_if_t<(is_pointer_t(typename MapU_t::mapped_type)), G4int> = 0>
  this_type& operator+=(const G4VTHitsMap<U, MapU_t>& right) const
  {
    MapU_t* aHitsMap = right.GetMap();
    for (auto itr = aHitsMap->begin(); itr != aHitsMap->end(); ++itr)
      add<U>(itr->first, *(itr->second));
    return static_cast<this_type&>(*(const_cast<this_type*>(this)));
  }
  //------------------------------------------------------------------------//
  template <typename U, typename MapU_t,
    enable_if_t<! (is_pointer_t(typename MapU_t::mapped_type)), G4int> = 0>
  this_type& operator+=(const G4VTHitsMap<U, MapU_t>& right) const
  {
    MapU_t* aHitsMap = right.GetMap();
    for (auto itr = aHitsMap->begin(); itr != aHitsMap->end(); ++itr)
      add<U>(itr->first, itr->second);
    return static_cast<this_type&>(*(const_cast<this_type*>(this)));
  }
  //------------------------------------------------------------------------//

 public:
  //------------------------------------------------------------------------//
  //  Insert a hit object. Total number of hit objects stored in this
  //  map is returned.
  //------------------------------------------------------------------------//
  //  Standard vector overload for any type
  //      assumes type T has overload of += operator for U
  //------------------------------------------------------------------------//
  template <typename U = T, enable_if_t<is_same_t(U, T), G4int> = 0>
  std::size_t add(const G4int& key, U*& aHit) const
  {
    vector_type* theHitsVector = GetVector(key);
    _add(theHitsVector, key, *aHit);
    return theHitsVector->size();
  }
  //------------------------------------------------------------------------//
  //  Overload for different types
  //------------------------------------------------------------------------//
  template <typename U = T, enable_if_t<! is_same_t(U, T), G4int> = 0>
  std::size_t add(const G4int& key, U*& aHit) const
  {
    vector_type* theHitsVector = GetVector(key);
    store_type hit = allocate();
    get_reference(hit) += *aHit;
    _add(theHitsVector, key, *hit);
    return theHitsVector->size();
  }
  //------------------------------------------------------------------------//
  //  Overload for same type
  //------------------------------------------------------------------------//
  template <typename U = T, enable_if_t<is_same_t(U, T), G4int> = 0>
  std::size_t add(const G4int& key, U& aHit) const
  {
    vector_type* theHitsVector = GetVector(key);
    _add(theHitsVector, key, aHit);
    return theHitsVector->size();
  }
  //------------------------------------------------------------------------//
  //  Overload for different types
  //------------------------------------------------------------------------//
  template <typename U = T, enable_if_t<! is_same_t(U, T), G4int> = 0>
  std::size_t add(const G4int& key, U& aHit) const
  {
    vector_type* theHitsVector = GetVector(key);
    _add(theHitsVector, key, aHit);
    return theHitsVector->size();
  }
  //------------------------------------------------------------------------//

 public:
  //------------------------------------------------------------------------//
  //  Set a hit object. Total number of hit objects stored in this
  //  map is returned.
  //------------------------------------------------------------------------//
  //  Overload for same type T
  //------------------------------------------------------------------------//
  template <typename U = T, enable_if_t<is_same_t(U, T), G4int> = 0>
  inline std::size_t set(const G4int& key, U*& aHit) const
  {
    vector_type* theHitsVector = GetVector(key);
    _assign(theHitsVector, key, aHit);
    return theHitsVector->size();
  }
  //------------------------------------------------------------------------//
  //  Overload for different types
  //------------------------------------------------------------------------//
  template <typename U = T, enable_if_t<! is_same_t(U, T), G4int> = 0>
  inline std::size_t set(const G4int& key, U*& aHit) const
  {
    vector_type* theHitsVector = GetVector(key);
    store_type hit = allocate();
    get_reference(hit) += *aHit;
    _assign(theHitsVector, key, hit);
    return theHitsVector->size();
  }
  //------------------------------------------------------------------------//

 public:
  //------------------------------------------------------------------------//
  //  Set a hit object. Total number of hit objects stored in this
  //  map is returned.
  //------------------------------------------------------------------------//
  //  Overload for same type T
  //------------------------------------------------------------------------//
  template <typename U = T, enable_if_t<is_same_t(U, T), G4int> = 0>
  inline std::size_t set(const G4int& key, U& aHit) const
  {
    vector_type* theHitsVector = GetVector(key);
    _assign(theHitsVector, key, &aHit);
    return theHitsVector->size();
  }
  //------------------------------------------------------------------------//
  //  Overload for different types
  //------------------------------------------------------------------------//
  template <typename U = T, enable_if_t<! is_same_t(U, T), G4int> = 0>
  inline std::size_t set(const G4int& key, U& aHit) const
  {
    vector_type* theHitsVector = GetVector(key);
    store_type hit = allocate();
    get_reference(hit) += aHit;
    _assign(theHitsVector, key, &aHit);
    return theHitsVector->size();
  }
  //------------------------------------------------------------------------//

 public:
  //------------------------------------------------------------------------//
  T* at(G4int key) const
  {
    vector_type* theHitsVector = GetVector();
    resize(theHitsVector, key);
    return &get_reference((*theHitsVector)[key]);
  }
  //------------------------------------------------------------------------//
  T* operator[](G4int key) const
  {
    vector_type* theHitsVector = GetVector();
    resize(theHitsVector, key);
    return &get_reference((*theHitsVector)[key]);
  }
  //------------------------------------------------------------------------//

 protected:
  template <typename U = store_type, enable_if_t<(is_pointer_t(U)), G4int> = 0>
  void resize(vector_type*& theHitsVector, const G4int& key) const
  {
    // ensure the proper size
    if (key >= (G4int)theHitsVector->size()) theHitsVector->resize(key + 1, null());

    // if null pointer for vector entry: allocate
    if (! theHitsVector->at(key)) {
      store_type init = allocate();
      _assign(theHitsVector, key, init);
    }
  }

  template <typename U = store_type, enable_if_t<(! is_pointer_t(U)), G4int> = 0>
  void resize(vector_type*& theHitsVector, const G4int& key) const
  {
    // ensure the proper size
    if (key >= (G4int)theHitsVector->size()) theHitsVector->resize(key + 1, null());
  }

  vector_type* GetVector(const G4int& key) const
  {
    vector_type* theHitsVector = GetVector();
    resize(theHitsVector, key);
    return theHitsVector;
  }

  //------------------------------------------------------------------------//
  //  Assign/Add when the storage type is pointer
  //      assumes type T has overload of += operator for U
  //------------------------------------------------------------------------//
  template <typename U = store_type, enable_if_t<is_pointer_t(U), G4int> = 0>
  void _assign(vector_type*& theHitsVector, const G4int& key, T& val) const
  {
    delete (*theHitsVector)[key];
    *(*theHitsVector)[key] = val;
  }

  template <typename U = store_type, enable_if_t<is_pointer_t(U), G4int> = 0>
  void _assign(vector_type*& theHitsVector, const G4int& key, T*& val) const
  {
    delete (*theHitsVector)[key];
    (*theHitsVector)[key] = val;
  }

  template <typename U = store_type, enable_if_t<is_pointer_t(U), G4int> = 0>
  void _add(vector_type*& theHitsVector, const G4int& key, T& val) const
  {
    *(*theHitsVector)[key] += val;
  }

  template <typename U = store_type, enable_if_t<is_pointer_t(U), G4int> = 0>
  void _add(vector_type*& theHitsVector, const G4int& key, T*& val) const
  {
    *(*theHitsVector)[key] += *val;
  }

  template <typename V, typename U = store_type, enable_if_t<is_pointer_t(U), G4int> = 0>
  void _add(vector_type*& theHitsVector, const G4int& key, V& val) const
  {
    *(*theHitsVector)[key] += val;
  }

  template <typename V, typename U = store_type, enable_if_t<is_pointer_t(U), G4int> = 0>
  void _add(vector_type*& theHitsVector, const G4int& key, V*& val) const
  {
    *(*theHitsVector)[key] += *val;
  }

  template <typename U = store_type, enable_if_t<is_pointer_t(U), G4int> = 0>
  T& get(U& val) const
  {
    return *val;
  }

  template <typename U = store_type, enable_if_t<is_pointer_t(U), G4int> = 0>
  void delete_contents(vector_type*& theHitsVector) const
  {
    for (auto itr = theHitsVector->begin(); itr != theHitsVector->end(); ++itr)
      delete *itr;
  }

  template <typename U = store_type, enable_if_t<is_pointer_t(U), G4int> = 0>
  T& get_reference(U& val) const
  {
    return *val;
  }

  //------------------------------------------------------------------------//
  //  Assign/Add when the storage type is pointer
  //      assumes type T has overload of += operator for U
  //------------------------------------------------------------------------//
  template <typename U = store_type, enable_if_t<! is_pointer_t(U), G4int> = 0>
  void _assign(vector_type*& theHitsVector, const G4int& key, T& val) const
  {
    (*theHitsVector)[key] = val;
  }

  template <typename U = store_type, enable_if_t<! is_pointer_t(U), G4int> = 0>
  void _assign(vector_type*& theHitsVector, const G4int& key, T*& val) const
  {
    delete (*theHitsVector)[key];
    (*theHitsVector)[key] = *val;
  }

  template <typename U = store_type, enable_if_t<! is_pointer_t(U), G4int> = 0>
  void _add(vector_type*& theHitsVector, const G4int& key, T& val) const
  {
    (*theHitsVector)[key] += val;
  }

  template <typename U = store_type, enable_if_t<! is_pointer_t(U), G4int> = 0>
  void _add(vector_type*& theHitsVector, const G4int& key, T*& val) const
  {
    (*theHitsVector)[key] += *val;
  }

  template <typename V, typename U = store_type, enable_if_t<! is_pointer_t(U), G4int> = 0>
  void _add(vector_type*& theHitsVector, const G4int& key, V& val) const
  {
    (*theHitsVector)[key] += val;
  }

  template <typename V, typename U = store_type, enable_if_t<! is_pointer_t(U), G4int> = 0>
  void _add(vector_type*& theHitsVector, const G4int& key, V*& val) const
  {
    (*theHitsVector)[key] += *val;
  }

  template <typename U = store_type, enable_if_t<! is_pointer_t(U), G4int> = 0>
  void delete_contents(vector_type*&) const
  {}

  template <typename U = store_type, enable_if_t<! is_pointer_t(U), G4int> = 0>
  T& get_reference(U& val) const
  {
    return val;
  }

#undef is_same_t
#undef is_fundamental_t
#undef is_std_map_t
#undef is_std_mmap_t
#undef is_std_uomap_t
#undef is_std_uommap_t
#undef is_map_t
#undef is_pointer_t
#undef scast
};

//============================================================================//

template <typename T, typename Vector_t>
G4VTHitsVector<T, Vector_t>::G4VTHitsVector(G4int init_size)
{
  theCollection = static_cast<void*>(new Vector_t);
  if (init_size > 0) {
    vector_type* theHitsVector = GetVector();
    resize(theHitsVector, init_size - 1);
  }
}

//============================================================================//

template <typename T, typename Vector_t>
G4VTHitsVector<T, Vector_t>::G4VTHitsVector(G4String detName, G4String colNam, G4int init_size)
  : G4HitsCollection(detName, colNam)
{
  theCollection = static_cast<void*>(new Vector_t);
  if (init_size > 0) {
    vector_type* theHitsVector = GetVector();
    resize(theHitsVector, init_size - 1);
  }
}

//============================================================================//

template <typename T, typename Vector_t>
G4VTHitsVector<T, Vector_t>::~G4VTHitsVector()
{
  vector_type* theHitsVector = GetVector();
  delete_contents(theHitsVector);
  delete theHitsVector;
}

//============================================================================//

template <typename T, typename Vector_t>
G4bool G4VTHitsVector<T, Vector_t>::operator==(const G4VTHitsVector<T, Vector_t>& right) const
{
  return (collectionName == right.collectionName);
}

//============================================================================//

template <typename T, typename Vector_t>
typename G4VTHitsVector<T, Vector_t>::map_t* G4VTHitsVector<T, Vector_t>::GetMap() const
{
  G4ThreadLocalStatic auto  theHitsMap = new map_t();
  theHitsMap->clear();
  vector_type* theHitsVector = GetVector();
  for (std::size_t i = 0; i < theHitsVector->size(); ++i) {
    store_type& _obj = (*theHitsVector)[i];
    if (_obj) (*theHitsMap)[i] = _obj;
  }
  return theHitsMap;
}
//============================================================================//

template <typename T, typename Vector_t>
void G4VTHitsVector<T, Vector_t>::DrawAllHits()
{
  ;
}

//============================================================================//

template <typename T, typename Vector_t>
void G4VTHitsVector<T, Vector_t>::PrintAllHits()
{
  G4cout << "G4THitsVector " << SDname << " / " << collectionName << " --- " << entries()
         << " entries" << G4endl;
  /*----- commented out for the use-case where <T> cannot be initialized
          to be zero or does not support += operator.
   Vector_t * theHitsVector = GetVector();
   typename Vector_t::iterator itr = theHitsVector->begin();
   T sum = 0.;
   for(; itr != theHitsVector->end(); ++itr) {
    G4cout << "  " << itr->first << " : "
           << *(itr->second) << G4endl;
    sum += *(itr->second);
   }
   G4cout << "             Total : " << sum << G4endl;
  ----------------------------------------------------------------------*/
}

//============================================================================//

template <typename T, typename Vector_t>
void G4VTHitsVector<T, Vector_t>::clear()
{
  vector_type* theHitsVector = GetVector();
  delete_contents(theHitsVector);
  theHitsVector->clear();
}

//============================================================================//
//                                                                            //
//                                                                            //
//                      Helpers for different map types                       //
//                                                                            //
//                                                                            //
//============================================================================//

template <typename _Tp>
class G4THitsVector : public G4VTHitsVector<_Tp, std::vector<_Tp*>>
{
 public:
  using parent_type = G4VTHitsVector<_Tp, std::vector<_Tp *>>;

 public:
  G4THitsVector(G4int init_size = 0) : parent_type(init_size) {}
  G4THitsVector(G4String detName, G4String colName, G4int init_size = 0)
    : parent_type(detName, colName, init_size)
  {}

  using parent_type::operator+=;
  using parent_type::operator==;
  using parent_type::operator[];
  using parent_type::add;
  using parent_type::begin;
  using parent_type::cbegin;
  using parent_type::cend;
  using parent_type::clear;
  using parent_type::DrawAllHits;
  using parent_type::end;
  using parent_type::entries;
  using parent_type::GetHit;
  using parent_type::GetMap;
  using parent_type::GetSize;
  using parent_type::GetVector;
  using parent_type::PrintAllHits;
  using parent_type::set;
};

//============================================================================//

template <typename _Tp>
class G4THitsDeque : public G4VTHitsVector<_Tp, std::deque<_Tp*>>
{
 public:
  using parent_type = G4VTHitsVector<_Tp, std::deque<_Tp*>>;

 public:
  G4THitsDeque(G4int init_size = 0) : parent_type(init_size) {}
  G4THitsDeque(G4String detName, G4String colName, G4int init_size = 0)
    : parent_type(detName, colName, init_size)
  {}

  using parent_type::operator+=;
  using parent_type::operator==;
  using parent_type::operator[];
  using parent_type::add;
  using parent_type::begin;
  using parent_type::cbegin;
  using parent_type::cend;
  using parent_type::clear;
  using parent_type::DrawAllHits;
  using parent_type::end;
  using parent_type::entries;
  using parent_type::GetHit;
  using parent_type::GetMap;
  using parent_type::GetSize;
  using parent_type::GetVector;
  using parent_type::PrintAllHits;
  using parent_type::set;
};

//============================================================================//

#endif
