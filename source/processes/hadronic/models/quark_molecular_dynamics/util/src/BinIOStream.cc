//
// $RCSfile: BinIOStream.cc,v $
//
// $Revision: 1.3 $
// $Date: 2003-06-19 14:43:45 $
// $Author: gunter $
// $Locker:  $
// $State: Exp $
// DOSfile: biostream.cpp
// UNIXfile: BinIOStream.C
//
// $Log: not supported by cvs2svn $
// Revision 1.2  1999/12/15 14:53:32  gunter
// Modifications to have a ISO C++ compliant source code
//
// Revision 1.1  1999/12/09 11:50:39  sscherer
// adding quark_molecular_dynamics
//
// Revision 1.1.1.1  1998/09/22 16:31:05  mhofmann
// U++ V0.9
//
// Revision 1.1.1.1  1997/07/17 12:54:53  mhofmann
// Initial import
//
// Revision 1.1.1.1  1996/10/04 14:37:49  mhofmann
// Phase Transition Project
//
//

static char RCS_Id[] = 
  "@(#) $Id: BinIOStream.cc,v 1.3 2003-06-19 14:43:45 gunter Exp $";

#include <iostream>
#include <stdlib.h>

#ifdef __MSDOS__
  #include "biostream.h"
#else
  #include "BinIOStream.hh"
#endif

#ifdef NO_INLINE
  #define inline
  #ifdef __MSDOS__
    #include "biostream.icp"
  #else
    #include "BinIOStream.icc"
  #endif
  #undef inline
#endif
     
static const char identification[] = "@(#)BinIOStream";

BinOStream::BinOStream(std::ostream& o, int) 
  : os(&o)
{
  for( int i = 0; i < 16; i++ ) {
    o << identification[i];
  }
}

BinIStream::BinIStream(std::istream& i, int) 
  : is(&i)
{
  char c; 
  for( int j = 0; j < 16; j++ ) {
    i >> c;
    if( c != identification[j] ) {
      exit(1);
    };
  }
}



