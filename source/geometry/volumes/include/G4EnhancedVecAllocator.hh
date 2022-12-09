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
// Class Description:
//
// A class for fast allocation of STL vectors through a static pool.
// It's meant to be used as alternative allocator for STL vectors.
       
//      ---------------- G4EnhancedVecAllocator ----------------
//
// Original author: X.Dong (NorthEastern Univ.), November 2009
// Reviewed implementation: G.Cosmo (CERN), December 2009
// ------------------------------------------------------------

#ifndef G4EnhancedVecAllocator_h
#define G4EnhancedVecAllocator_h 1

#include "G4Types.hh"

typedef struct
{
  G4int isAllocated;
  char *address;
} G4ChunkType;

typedef struct
{
  std::size_t size;
  G4int totalspace;
  G4ChunkType *preAllocated;
} G4ChunkIndexType;

class G4AllocStats
{
  // --------------------------------------------------------------------
  // Utility class, placeholder for global data on allocation.
  // Initialisation to zero of the data below *must* be added ONCE only
  // directly in the client code, where this allocator is to be applied
  // --------------------------------------------------------------------

  public:

    static G4ThreadLocal G4ChunkIndexType * allocStat;
    static G4ThreadLocal G4int totSpace;
    static G4ThreadLocal G4int numCat;
};

template<typename _Tp>
class G4EnhancedVecAllocator : public std::allocator<_Tp>
{
  public:

    template<typename _Tp1>
    struct rebind { typedef G4EnhancedVecAllocator<_Tp1> other; };

    G4EnhancedVecAllocator() {;}

    G4EnhancedVecAllocator(const G4EnhancedVecAllocator<_Tp>&)
      : std::allocator<_Tp>() {;}

    template<typename _Tp1>
    G4EnhancedVecAllocator(const G4EnhancedVecAllocator<_Tp1>&)
      : std::allocator<_Tp>() {;}

    ~G4EnhancedVecAllocator() {;}

    // override allocate / deallocate
    //
    void deallocate(_Tp* _Ptr, std::size_t _Count);
#ifdef __IBMCPP__
    _Tp* allocate(std::size_t _Count, void * const hint = 0);  // IBM AIX
#else
    _Tp* allocate(std::size_t _Count);
#endif
};

// ------------------------------------------------------------
// Inline implementations
// ------------------------------------------------------------

// ************************************************************
// deallocate
// ************************************************************
//
template<typename _Tp>
void G4EnhancedVecAllocator<_Tp>::deallocate(_Tp* _Ptr, std::size_t _Count)
{
  G4int found = -1;
  if (G4AllocStats::allocStat != 0)
  {
    for (auto j = 0 ; j < G4AllocStats::numCat ; ++j)
    {
      if (G4AllocStats::allocStat[j].size == (_Count * sizeof(_Tp)))
      {
        found = j;
        break;
      }
    }
  }
  // assert(found != -1);

  G4ChunkIndexType& chunk = G4AllocStats::allocStat[found];

  for (auto k = 0; k < chunk.totalspace; ++k)
  {
    if ( (chunk.preAllocated[k]).address == ((char *) _Ptr))
    {
      // assert((chunk.preAllocated[k]).isAllocated==1);
      (chunk.preAllocated[k]).isAllocated = 0;
      return;
    }
  }
}

// ************************************************************
// allocate
// ************************************************************
//
#ifdef __IBMCPP__
template<typename _Tp>
_Tp* G4EnhancedVecAllocator<_Tp>::allocate(std::size_t _Count, void * const hint)
#else
template<typename _Tp>
_Tp* G4EnhancedVecAllocator<_Tp>::allocate(std::size_t _Count)
#endif
{
  std::size_t totalsize = _Count * sizeof(_Tp);

  G4int found = -1;
  if (G4AllocStats::allocStat != 0)
  {
    for (auto j = 0 ; j < G4AllocStats::numCat ; ++j)
    {
      if (G4AllocStats::allocStat[j].size == totalsize)
      {
        found = j;
        break;
      } 
    }
  }   

  if (found == -1)  // Find the new size
  {
    ++G4AllocStats::numCat;
    if (G4AllocStats::numCat > G4AllocStats::totSpace)
    {
      G4AllocStats::totSpace = G4AllocStats::totSpace + 128;
        // heuristic parameter for different sizes

      G4AllocStats::allocStat =
           (G4ChunkIndexType *) realloc(G4AllocStats::allocStat,
           sizeof(G4ChunkIndexType) * G4AllocStats::totSpace);
        // This value must be different than zero; otherwise means
        // failure in allocating extra space !
      // assert(G4AllocStats::allocStat != 0);
    }

    G4AllocStats::allocStat[G4AllocStats::numCat-1].size = totalsize;
    G4AllocStats::allocStat[G4AllocStats::numCat-1].totalspace = 0;
    G4AllocStats::allocStat[G4AllocStats::numCat-1].preAllocated = 0;

    found = G4AllocStats::numCat - 1;
    G4ChunkIndexType& chunk = G4AllocStats::allocStat[found];

    chunk.totalspace = 512;
      // heuristic for the number of STL vector instances

    chunk.preAllocated = (G4ChunkType *) realloc(chunk.preAllocated,
                         sizeof(G4ChunkType) * chunk.totalspace);
      // This value must be different than zero; otherwise means
      // failure in allocating extra space for pointers !
    // assert(chunk.preAllocated != 0);

    char *newSpace1 = (char *) malloc(totalsize * 512);
      // This pointer must be different than zero; otherwise means
      // failure in allocating extra space for instances !
    // assert(newSpace1 != 0);

    for (auto k = 0; k < 512 ; ++k)
    {
      (chunk.preAllocated[k]).isAllocated = 0;
      (chunk.preAllocated[k]).address = newSpace1+totalsize*k;
    }

    (chunk.preAllocated[0]).isAllocated = 1;
    return (_Tp*)((chunk.preAllocated[0]).address);
  }

  G4ChunkIndexType& chunk = G4AllocStats::allocStat[found];

  // assert(chunk.size == totalsize);

  for (auto k = 0; k < chunk.totalspace; ++k)
  {
    if ((chunk.preAllocated[k]).isAllocated == 0)
    { 
      (chunk.preAllocated[k]).isAllocated = 1;
      return (_Tp*)((chunk.preAllocated[k]).address);
    }
  }

  G4int originalchunknumber = chunk.totalspace;
      
  chunk.totalspace += 512;  // heuristic for the number of STL vector instances

  chunk.preAllocated = (G4ChunkType *) realloc(chunk.preAllocated,
                       sizeof(G4ChunkType) * chunk.totalspace);
    // This value must be different than zero; otherwise means
    // failure in allocating extra space for pointers !
  // assert(chunk.preAllocated != 0);

  char *newSpace = (char *) malloc(totalsize * 512);
    // This pointer must be different than zero; otherwise means
    // failure in allocating extra space for instances !
  // assert(newSpace != 0);

  for (auto k = 0; k < 512 ; ++k)
  {
    (chunk.preAllocated[originalchunknumber+k]).isAllocated = 0;
    (chunk.preAllocated[originalchunknumber+k]).address = newSpace+totalsize*k;
  }

  (chunk.preAllocated[originalchunknumber]).isAllocated = 1;

  return (_Tp*)((chunk.preAllocated[originalchunknumber]).address);
}

// ************************************************************
// operator==
// ************************************************************
//
template<typename _T1, typename _T2>
inline G4bool operator==(const G4EnhancedVecAllocator<_T1>&,
                         const G4EnhancedVecAllocator<_T2>&)
{ return true; }

// ************************************************************
// operator!=
// ************************************************************
//
template<typename _T1, typename _T2>
inline G4bool operator!=(const G4EnhancedVecAllocator<_T1>&,
                         const G4EnhancedVecAllocator<_T2>&)
{ return false; }

#endif
