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
//---------------------------------------------------------------
//
// G4DNAModelSubType.hh
//
// Class Description:
//   This is an enumerator to define sub-type of DNA models
//
// Creation date: 18.05.2018
// Modifications:
//
//---------------------------------------------------------------

#ifndef G4DNAModelSubType_h
#define G4DNAModelSubType_h 1

enum G4DNAModelSubType 
{
  fDNAUnknownModel = 0,
  fRitchie1994eSolvation = 1, 
  fTerrisol1990eSolvation = 2, 
  fMeesungnoen2002eSolvation = 3,
  fKreipl2009eSolvation = 4,
  fMeesungnoensolid2002eSolvation = 5
};

#endif
