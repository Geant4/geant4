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
// * GEANT4 collaboration and of QinetiQ Ltd,  subject DEFCON 705 IPR *
// * conditions.                                                      *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4VFacet.cc,v 1.2 2006-01-30 14:39:53 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4VFacet.hh
//
// Date:                15/06/2005
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:            C/MAT/N03517
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "G4VFacet.hh"
#include "globals.hh"

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
