// $Id: G4VPHitsCollectionIO.cc,v 1.2 2002-12-04 10:25:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VPHitsCollectionIO.cc
//
// History:
//   '01.08.16  Youhei Morita  Initial creation

#include "G4VPHitsCollectionIO.hh"

// Implementation of Constructor #1
G4VPHitsCollectionIO::G4VPHitsCollectionIO( G4std::string detName,
                                            G4std::string colName )
 : m_verbose(0), f_detName(detName), f_colName(colName)
{}

// Implementation of operator== 
G4bool G4VPHitsCollectionIO::operator== (const G4VPHitsCollectionIO& right) const
{
  return ( (f_detName == right.f_detName) &&
           (f_colName == right.f_colName) );
}

// End of G4VPHitsCollectionIO.cc

