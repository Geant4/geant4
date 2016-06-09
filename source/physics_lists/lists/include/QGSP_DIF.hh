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
// $Id: QGSP_DIF.hh,v 1.1 2007/11/13 10:20:32 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   QGSP_DIF
//
// Author: 2007  G.Folger
//           created from QGSP  originally by J.P. Wellisch
//
// Modified:
//
//----------------------------------------------------------------------------
#ifndef TQGSP_DIF_h
#define TQGSP_DIF_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TQGSP_DIF: public T
{
public:
  TQGSP_DIF(G4int ver = 1);
  virtual ~TQGSP_DIF();
  
public:
  // SetCuts() 
  virtual void SetCuts();
  
private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };

};
#include "QGSP_DIF.icc"
typedef TQGSP_DIF<G4VModularPhysicsList> QGSP_DIF;

// 2002 by J.P. Wellisch

#endif



