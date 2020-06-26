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
// ---------------------------------------------------------------
//  Tasking class implementation
//
// TaskAllocatorPool
//
// Implementation file
//
// Author: G.Cosmo, November 2000
//

#include "PTL/TaskAllocatorPool.hh"

using namespace PTL;

// ************************************************************
// TaskAllocatorPool constructor
// ************************************************************
//
TaskAllocatorPool::TaskAllocatorPool(unsigned int sz)
: esize((sz < sizeof(PoolLink)) ? sizeof(PoolLink) : sz)
, csize((sz < 1024 / 2 - 16) ? (1024 - 16) : (sz * 10 - 16))
, chunks(nullptr)
, head(nullptr)
, nchunks(0)
{}

// ************************************************************
// TaskAllocatorPool copy constructor
// ************************************************************
//
TaskAllocatorPool::TaskAllocatorPool(const TaskAllocatorPool& right)
: esize(right.esize)
, csize(right.csize)
, chunks(right.chunks)
, head(right.head)
, nchunks(right.nchunks)
{}

// ************************************************************
// TaskAllocatorPool operator=
// ************************************************************
//
TaskAllocatorPool&
TaskAllocatorPool::operator=(const TaskAllocatorPool& right)
{
    if(&right == this)
        return *this;
    chunks  = right.chunks;
    head    = right.head;
    nchunks = right.nchunks;
    return *this;
}

// ************************************************************
// TaskAllocatorPool destructor
// ************************************************************
//
TaskAllocatorPool::~TaskAllocatorPool() { Reset(); }

// ************************************************************
// Reset
// ************************************************************
//
void
TaskAllocatorPool::Reset()
{
    // Free all chunks
    //
    PoolChunk* n = chunks;
    PoolChunk* p = nullptr;
    while(n)
    {
        p = n;
        n = n->next;
        delete p;
    }
    head    = nullptr;
    chunks  = nullptr;
    nchunks = 0;
}

// ************************************************************
// Grow
// ************************************************************
//
void
TaskAllocatorPool::Grow()
{
    // Allocate new chunk, organize it as a linked list of
    // elements of size 'esize'
    //
    PoolChunk* n = new PoolChunk(csize);
    n->next      = chunks;
    chunks       = n;
    nchunks++;

    const int nelem = csize / esize;
    char*     start = n->mem;
    char*     last  = &start[(nelem - 1) * esize];
    for(char* p = start; p < last; p += esize)
    {
        reinterpret_cast<PoolLink*>(p)->next = reinterpret_cast<PoolLink*>(p + esize);
    }
    reinterpret_cast<PoolLink*>(last)->next = nullptr;
    head                                    = reinterpret_cast<PoolLink*>(start);
}
