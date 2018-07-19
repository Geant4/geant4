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
// $Id: G4NuclearFormfactorType.hh 84082 2014-10-06 14:33:09Z urban $
//
//---------------------------------------------------------------
//
// G4NuclearFormfactorType.hh
//
// Class Description:
//   This is an enumerator to define type nuclear formfactor
//   parameterisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 9.03.2016
// Modifications:
//
//---------------------------------------------------------------

#ifndef G4NuclearFormfactorType_h
#define G4NuclearFormfactorType_h 1

enum G4NuclearFormfactorType
{
  fNoneNF = 0,
  fExponentialNF,
  fGaussianNF,
  fFlatNF
};
#endif


