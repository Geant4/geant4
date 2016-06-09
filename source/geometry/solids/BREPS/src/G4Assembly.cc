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
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4Assembly.cc,v 1.7 2010-07-07 14:45:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Assembly.cc
//
// ----------------------------------------------------------------------

#include "G4Assembly.hh"

G4Assembly::G4Assembly()
  : numberOfSolids(0)
{
  //  ReadSTEPFile();
  //  CopySTEPData();  
}

G4Assembly::~G4Assembly()
{
  G4PlacedSolid* a = 0;
  
  // Remove placedVec and delete all its contents
  while (placedVec.size()>0)
  {
    a = placedVec.back();
    placedVec.pop_back();
    for (G4PlacedVector::iterator i=placedVec.begin(); i!=placedVec.end();)
    {
      if (*i==a)
      {
	i = placedVec.erase(i);
      }
      else
      {
	++i;
      }
    } 
    if ( a )  { delete a; }
  } 
}

void G4Assembly::SetPlacedVector(G4PlacedVector& pVec)
{
  numberOfSolids = pVec.size();
  
  for(G4int a=0;a<numberOfSolids;a++)
  {
    placedVec.push_back( pVec[a]);
  }
}
