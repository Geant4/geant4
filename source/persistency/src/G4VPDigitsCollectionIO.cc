// $Id: G4VPDigitsCollectionIO.cc,v 1.1 2002-11-24 13:45:25 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VPDigitsCollectionIO.cc
//
// History:
//   '01.08.16  Youhei Morita  Initial creation

#include "G4VPDigitsCollectionIO.hh"

// Implementation of Constructor #1
G4VPDigitsCollectionIO::G4VPDigitsCollectionIO( std::string detName, std::string colName )
 : m_verbose(0), f_detName(detName), f_colName(colName)
{}

// Implementation of operator== 
bool G4VPDigitsCollectionIO::operator== (const G4VPDigitsCollectionIO& right) const
{
  return ( (f_detName == right.f_detName) &&
           (f_colName == right.f_colName) );
}

// End of G4VPDigitsCollectionIO.cc

