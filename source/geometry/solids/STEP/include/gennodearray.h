

//



//
// $Id: gennodearray.h,v 1.2 1999-05-21 20:20:40 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef gennodearray_h
#define gennodearray_h

/*
* NIST Utils Class Library
* clutils/gennode.h
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 

/*
 * GenNodeArray - dynamic array object of GenericNodes.
 * the array part of this is copied from Unidraws UArray - dynamic array object
 * Copyright (c) 1990 Stanford University
 */

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

#include <string.h>
#include <stdlib.h> // to get bcopy for CenterLine

#include <gennode.h>

	// the initial size of the array
#define ARRAY_DEFAULT_SIZE (1024)

//////////////////////////////////////////////////////////////////////////////
// GenNodeArray
// If you delete this object, it does not delete the entries it points to.
// If you want it to delete the entries it points to you need to call
// DeleteEntries().
//////////////////////////////////////////////////////////////////////////////

class GenNodeArray 
{
public:
    GenNodeArray(int defaultSize = ARRAY_DEFAULT_SIZE);
    virtual ~GenNodeArray();

    GenericNode*& operator[](int index);
    virtual int Index(GenericNode* gn);
    virtual int Index(GenericNode** gn);

    int Count();

    virtual void Append(GenericNode* gn);
    virtual int Insert(GenericNode* gn);
    virtual int Insert(GenericNode* gn, int index);
    virtual void Remove(int index);
    virtual void ClearEntries();
    virtual void DeleteEntries();

protected:
    virtual void Check(int index);

    GenericNode** _buf;	// the array
    int _bufsize;	// the possible number of entries in the array
    int _count;		// the number of entries in the array
};

//////////////////////////////////////////////////////////////////////////////
// class GenNodeArray inline public functions
//////////////////////////////////////////////////////////////////////////////

inline GenNodeArray::GenNodeArray (int defaultSize)
{
    _bufsize = defaultSize;
    _buf = new GenericNode*[_bufsize];
    memset(_buf, 0, _bufsize*sizeof(GenericNode*));
    _count = 0;
}

inline GenNodeArray::~GenNodeArray ()
{

//    int i;
	// this is dangerous because several things point at these nodes
	// also whatever is derived from this thing might do this
//    for(i = 0; i < _count; i++)
//	delete _buf[i];
    delete [] _buf;
}

inline GenericNode*& GenNodeArray::operator[] (int index) 
{
    Check(index);
    return _buf[index];
}

inline int GenNodeArray::Index (GenericNode** gn)
{
    return ((gn - _buf) / sizeof(GenericNode*));
}

inline void GenNodeArray::Append(GenericNode* gn)
{    
    Insert(gn, _count); 
}

inline int GenNodeArray::Insert(GenericNode* gn)
{
    return Insert(gn, _count); 
}

inline int GenNodeArray::Count ()
{
    return _count;
}

#endif
