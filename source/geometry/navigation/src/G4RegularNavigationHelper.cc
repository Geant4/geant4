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
// $Id: G4RegularNavigationHelper.cc,v 1.1 2009-01-27 09:31:29 gcosmo Exp $
// GEANT4 tag $ Name:$
//
// class G4RegularNavigationHelper implementation
//
// Author: Pedro Arce, November 2008
//
// --------------------------------------------------------------------

#include "G4RegularNavigationHelper.hh"

std::vector< std::pair<G4int,G4double> >
  G4RegularNavigationHelper::theStepLengths;

// --------------------------------------------------------------------
//
G4RegularNavigationHelper::G4RegularNavigationHelper()
{
}

// --------------------------------------------------------------------
//
G4RegularNavigationHelper::~G4RegularNavigationHelper()
{
}

// --------------------------------------------------------------------
//
void G4RegularNavigationHelper::ClearStepLengths()
{
  theStepLengths.clear();
}

// --------------------------------------------------------------------
//
void G4RegularNavigationHelper::AddStepLength( G4int copyNo, G4double slen )
{
  theStepLengths.push_back( std::pair<G4int,G4double>(copyNo,slen) );
}
