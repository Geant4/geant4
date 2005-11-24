// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		G4VFacet.hh
//
// Date:		15/06/2005
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:		C/MAT/N03517
//
// This software is the intelectual property of QinetiQ Ltd, subject
// DEFCON 705 IPR conditions.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DISCLAIMER
// ----------
//
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
//
//
///////////////////////////////////////////////////////////////////////////////
//
//
#include "G4VFacet.hh"
#include "globals.hh"
using namespace std;
///////////////////////////////////////////////////////////////////////////////
//
std::ostream &G4VFacet::StreamInfo(std::ostream &os) const
{
  os <<G4endl;
  os <<"***********************************************************************"
     <<G4endl;
  os <<"FACET TYPE       = " <<geometryType <<G4endl;
  os <<"ABSOLUTE VECTORS = " <<G4endl;
  os <<"P0               = " <<P0 <<G4endl;
  for (G4ThreeVectorList::const_iterator it=P.begin(); it!=P.end(); it++)
    os <<"P[" <<it-P.begin()+1 <<"]      = " <<*it <<G4endl;
    
  os <<"RELATIVE VECTORS = " <<G4endl;
  for (G4ThreeVectorList::const_iterator it=E.begin(); it!=E.end(); it++)
    os <<"E[" <<it-E.begin()+1 <<"]      = " <<*it <<G4endl;
  
  os <<"***********************************************************************"
     <<G4endl;
  
  return os;
}
///////////////////////////////////////////////////////////////////////////////
//
