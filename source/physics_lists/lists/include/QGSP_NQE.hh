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
// $Id: QGSP_NQE.hh,v 1.1 2007/04/26 14:47:11 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   QGSP_NQE
//
// Author: 2007 G.Folger
//     created from QGSP 
//
// Modified:
//
//----------------------------------------------------------------------------
#ifndef TQGSP_NQE_h
#define TQGSP_NQE_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TQGSP_NQE: public T
{
public:
  TQGSP_NQE(G4int ver = 1);
  virtual ~TQGSP_NQE();
  
public:
  // SetCuts() 
  virtual void SetCuts();
  
private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };

};
#include "QGSP_NQE.icc"
typedef TQGSP_NQE<G4VModularPhysicsList> QGSP_NQE;

// 2002 by J.P. Wellisch

#endif



