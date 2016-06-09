//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject DEFCON 705 IPR conditions.                               *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4VFacet.cc,v 1.3 2006/06/29 18:49:31 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
