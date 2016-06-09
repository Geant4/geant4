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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName: LHEP_EMV
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.12.2005 G.Folger: create from HadronPhysicsLHEP_GN
//
//----------------------------------------------------------------------------
//
#ifndef TLHEP_EMV_h
#define TLHEP_EMV_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VModularPhysicsList.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TLHEP_EMV: public T
{
public:
  TLHEP_EMV(G4int ver = 1);
  virtual ~TLHEP_EMV();
  
public:
  // SetCuts() 
  virtual void SetCuts();
private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };

};
#include "LHEP_EMV.icc"
typedef TLHEP_EMV<G4VModularPhysicsList> LHEP_EMV;

// 2002 by J.P. Wellisch

#endif



