#ifndef CSVofstream_h
#define CSVofstream_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		CSVofstream.hh
//
// Version:		0.a1
// Date:		06/02/99
// Author:		P R Truscott
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//
// DESCRIPTION
// -----------
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
////////////////////////////////////////////////////////////////////////////////
//
#include "fstream.h"

////////////////////////////////////////////////////////////////////////////////
//
class CSVofstream : public ofstream
{
  public:
    CSVofstream () : ofstream() {};
  //    CSVofstream (const char *name, int mode=std::ios::out, int prot=0664)
  //      : ofstream(name, mode, prot) {};
  CSVofstream (const char *name)
        : ofstream(name) {};
    ~CSVofstream () {};
};
////////////////////////////////////////////////////////////////////////////////
#endif
