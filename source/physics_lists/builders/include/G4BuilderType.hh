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
// $Id: G4BuilderType.hh 66241 2012-12-13 18:34:42Z gunter $
//
//---------------------------------------------------------------
//
// G4BuilderType.hh
//
// Class Description:
//   This is an enumerator to define sub-type of electro-magnetic
//   processes
//
// Creation date: 23.09.2008
// Modifications:
//
//---------------------------------------------------------------

#ifndef G4BuilderType_h
#define G4BuilderType_h 1

enum G4BuilderType 
{ 
  bUnknown = 0, 
  bTransportation = 1, 
  bElectromagnetic = 2, 
  bEmExtra = 3,
  bDecay = 4,
  bHadronElastic = 5, 
  bHadronInelastic = 6,
  bStopping = 7,
  bIons = 8
};

#endif
