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
// $Id: QBBC.hh 66892 2013-01-17 10:57:59Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:  QBBC
//
// Author: 11 April 2006 V. Ivanchenko
//
// Modified:
// 15.04.2007 set glauber=true (V.Ivanchenko)
//----------------------------------------------------------------------------
//
#ifndef QBBC_h
#define QBBC_h 1

#include "globals.hh"
#include "G4VModularPhysicsList.hh"

class QBBC : public G4VModularPhysicsList
{
public:

  QBBC(G4int ver = 1, const G4String& type = "QBBC");

  virtual ~QBBC();

  virtual void SetCuts();

private:

  // copy constructor and hide assignment operator
  QBBC(QBBC &);
  QBBC & operator=(const QBBC &right);

};

#endif



