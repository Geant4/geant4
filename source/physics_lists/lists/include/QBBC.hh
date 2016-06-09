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
// $Id: QBBC.hh,v 1.2 2007/04/16 11:57:40 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
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
#ifndef TQBBC_h
#define TQBBC_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TQBBC: public T
{
public:

  TQBBC(G4int ver = 1, const G4String& type = "QBBC", G4bool glauber = true);

  virtual ~TQBBC();
  
public:

  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};

#include "QBBC.icc"
typedef TQBBC<G4VModularPhysicsList> QBBC;

#endif



