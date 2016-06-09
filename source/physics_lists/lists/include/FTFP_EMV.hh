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
// $Id: FTFP_EMV.hh,v 1.1 2007/04/26 14:47:10 gunter Exp $
// GEANT4 tag $Name: geant4-08-03 $
//
//---------------------------------------------------------------------------
//
// ClassName: FTFP_EMV
//
// Author: 2007 G.Folger
//    Created as copy from FTFP
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef TFTFP_EMV_h
#define TFTFP_EMV_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TFTFP_EMV: public T
{
public:
  TFTFP_EMV(G4int ver = 1);
  virtual ~TFTFP_EMV();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};
#include "FTFP_EMV.icc"
typedef TFTFP_EMV<G4VModularPhysicsList> FTFP_EMV;

// 2002 by J.P. Wellisch

#endif



