//
// MIT License
// Copyright (c) 2020 Jonathan R. Madsen
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED
// "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
// LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
// -------------------------------------------------------------------
//      Tasking class header file
//
// Class description:
//
// Class implementing a memory pool for fast allocation and deallocation
// of memory chunks.  The size of the chunks for small allocated objects
// is fixed to 1Kb and takes into account of memory alignment; for large
// objects it is set to 10 times the object's size.
// The implementation is derived from: B.Stroustrup, The C++ Programming
// Language, Third Edition.

//           -------------- TaskAllocatorPool ----------------
//
// Author: G.Cosmo (CERN), November 2000
// -------------------------------------------------------------------

#pragma once

namespace PTL
{
class TaskAllocatorPool
{
public:
    explicit TaskAllocatorPool(unsigned int n = 0);
    // Create a pool of elements of size n
    ~TaskAllocatorPool();
    // Destructor. Return storage to the free store

    inline void* Alloc();
    // Allocate one element
    inline void Free(void* b);
    // Return an element back to the pool

    inline unsigned int Size() const;
    // Return storage size
    void Reset();
    // Return storage to the free store

    inline int GetNoPages() const;
    // Return the total number of allocated pages
    inline unsigned int GetPageSize() const;
    // Accessor for default page size
    inline void GrowPageSize(unsigned int factor);
    // Increase default page size by a given factor

private:
    TaskAllocatorPool(const TaskAllocatorPool& right);
    // Provate copy constructor
    TaskAllocatorPool& operator=(const TaskAllocatorPool& right);
    // Private equality operator

    struct PoolLink
    {
        PoolLink* next;
    };

    class PoolChunk
    {
    public:
        explicit PoolChunk(unsigned int sz)
        : size(sz)
        , mem(new char[size])
        , next(0)
        {
            ;
        }
        ~PoolChunk() { delete[] mem; }
        const unsigned int size;
        char*              mem;
        PoolChunk*         next;
    };

    void Grow();
    // Make pool larger

private:
    const unsigned int esize;
    unsigned int       csize;
    PoolChunk*         chunks;
    PoolLink*          head;
    int                nchunks;
};

//--------------------------------------------------------------------------------------//
// Inline implementation

// ************************************************************
// Alloc
// ************************************************************
//
inline void*
TaskAllocatorPool::Alloc()
{
    if(head == nullptr)
        Grow();
    PoolLink* p = head;  // return first element
    head        = p->next;
    return p;
}

// ************************************************************
// Free
// ************************************************************
//
inline void
TaskAllocatorPool::Free(void* b)
{
    PoolLink* p = static_cast<PoolLink*>(b);
    p->next     = head;  // put b back as first element
    head        = p;
}

// ************************************************************
// Size
// ************************************************************
//
inline unsigned int
TaskAllocatorPool::Size() const
{
    return nchunks * csize;
}

// ************************************************************
// GetNoPages
// ************************************************************
//
inline int
TaskAllocatorPool::GetNoPages() const
{
    return nchunks;
}

// ************************************************************
// GetPageSize
// ************************************************************
//
inline unsigned int
TaskAllocatorPool::GetPageSize() const
{
    return csize;
}

// ************************************************************
// GrowPageSize
// ************************************************************
//
inline void
TaskAllocatorPool::GrowPageSize(unsigned int sz)
{
    csize = (sz) ? sz * csize : csize;
}

}  // namespace PTL
