

//



//
// $Id: gennodearray.cc,v 1.2 1999-05-21 20:21:09 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST Utils Class Library
* clutils/gennodearray.h
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 

#include <gennode.h>
#include <gennodelist.h>
#include <gennodearray.h>

void 
GenNodeArray::Check (int index)
{
    GenericNode** newbuf;

    if (index >= _bufsize) {
	int oldBufSize = _bufsize;
        _bufsize = (index+1) * 2;
        newbuf = new GenericNode*[_bufsize];
	memset(newbuf, 0, _bufsize);
//	memset(newbuf[oldBufSize], 0, 
//		(_bufsize - oldBufSize)*sizeof(GenericNode*) );
//        bcopy(_buf, newbuf, _count*sizeof(GenericNode*));
// Josh L, 5/2/95
//        memcpy(newbuf, _buf, _count*sizeof(GenericNode*));
// Dave memcpy is not working since memory areas overlap
        memmove(newbuf, _buf, _count*sizeof(GenericNode*));
	delete [] _buf;
        _buf = newbuf;
    }
}

int 
GenNodeArray::Insert (GenericNode* gn, int index) 
{
    const GenericNode** spot;
    index = (index < 0) ? _count : index;

    if (index < _count) {
        Check(_count+1);
        spot = (const GenericNode**)&_buf[index];
//        bcopy(spot, spot+1, (_count - index)*sizeof(GenericNode*));
// Josh L, 5/2/95
//        memcpy(spot+1, spot, (_count - index)*sizeof(GenericNode*));
// Dave memcpy is not working since memory areas overlap
        memmove(spot+1, spot, (_count - index)*sizeof(GenericNode*));

    } else {
        Check(index);
        spot = (const GenericNode**)&_buf[index];
    }
    *spot = gn;
    ++_count;
    return index;
}

void 
GenNodeArray::Remove (int index) 
{
    if (0 <= index && index < _count) {
        --_count;
        const GenericNode** spot = (const GenericNode**)&_buf[index];
//        bcopy(spot+1, spot, (_count - index)*sizeof(GenericNode*));
// Josh L, 5/2/95
//        memcpy(spot, spot+1, (_count - index)*sizeof(GenericNode*));
// Dave memcpy is not working since memory areas overlap
        memmove(spot, spot+1, (_count - index)*sizeof(GenericNode*));
	_buf[_count] = 0;
    }
}

void GenNodeArray::ClearEntries ()
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "GenNodeArray::Clear()\n";
    int i;
    for(i = 0 ; i < _count; i++)
	_buf[i] = 0;
    _count = 0;
}

void GenNodeArray::DeleteEntries()
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "GenNodeArray::DeleteEntries()\n";
    int i;
    for(i = 0 ; i < _count; i++)
	delete (_buf[i]);
    _count = 0;
}


int GenNodeArray::Index (GenericNode* gn)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "GenNodeArray::Index()\n";
    for (int i = 0; i < _count; ++i) {
        if (_buf[i] == gn) {
            return i;
        }
    }
    return -1;
}

