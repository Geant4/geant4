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
// $Id: G4AllocatorPool.hh,v 1.2 2004/11/12 16:25:34 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// 
// -------------------------------------------------------------------
//      GEANT 4 class header file 
//
// Class description:
//
// Class implementing a memory pool for fast allocation and deallocation
// of memory chunks.  The size of the chunks for small allocated objects
// is fixed to 1Kb and takes into account of memory alignment; for large
// objects it is set to 10 times the object's size.
// The implementation is derived from: B.Stroustrup, The C++ Programming
// Language, Third Edition.

//           -------------- G4AllocatorPool ----------------
//
// Author: G.Cosmo (CERN), November 2000
// -------------------------------------------------------------------

#ifndef G4AllocatorPool_h
#define G4AllocatorPool_h 1

class G4AllocatorPool
{
  public:

    G4AllocatorPool( unsigned int n );
      // Create a pool of elements of size n
    ~G4AllocatorPool();
      // Destructor. Return storage to the free store

    G4AllocatorPool(const G4AllocatorPool& right);
      // Copy constructor

    inline void* Alloc();
      // Allocate one element
    inline void  Free( void* b );
      // Return an element back to the pool

    inline unsigned int  Size() const;
      // Return storage size
    void  Reset();
      // Return storage to the free store

  private:

    G4AllocatorPool& operator= (const G4AllocatorPool& right);
      // Private equality operator

    struct G4PoolLink
    {
      G4PoolLink* next;
    };
    class G4PoolChunk
    {
      public:
        G4PoolChunk(unsigned int sz)
          : size(sz), mem(new char[size]), next(0) {;}
        ~G4PoolChunk() { delete [] mem; }
        const unsigned int size;
        char* mem;
        G4PoolChunk* next;
    };

    void Grow();
      // Make pool larger

  private:

    G4PoolChunk* chunks;
    const unsigned int esize;
    const unsigned int csize;
    G4PoolLink* head;
    int nchunks;
};

// ------------------------------------------------------------
// Inline implementation
// ------------------------------------------------------------

// ************************************************************
// Alloc
// ************************************************************
//
inline void*
G4AllocatorPool::Alloc()
{
  if (head==0) Grow();
  G4PoolLink* p = head;  // return first element
  head = p->next;
  return p;
}

// ************************************************************
// Free
// ************************************************************
//
inline void
G4AllocatorPool::Free( void* b )
{
  G4PoolLink* p = static_cast<G4PoolLink*>(b);
  p->next = head;        // put b back as first element
  head = p;
}

// ************************************************************
// Size
// ************************************************************
//
inline unsigned int
G4AllocatorPool::Size() const
{
  return nchunks*csize;
}

#endif
