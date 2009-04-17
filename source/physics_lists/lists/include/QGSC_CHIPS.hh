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
// $Id: QGSC_CHIPS.hh,v 1.3 2009-04-17 15:29:19 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName: QGSC_CHIPS
//
// Author: 2009 M.Kosov (Development of the QGSC_QGSC physics list)
//
// Modified:
//
//-------------------------------------------------------------------------------
//  Short description: In tis file the aim is to get rid of LHEP, use the
//  QGSC model with the Energy Flow interface for CHIPS from E=0 for all
//  hadrons and gradually substitute QGSC by the native CHIPS (G4QCollision)
//  process starting from low energies (first nucleons, then mesons, after that
//  hyperons and antibaryons). The present "QGSC" model at low energies is
//  jus a temporary (far from ideal) interface to the CHIPS model. The nuclear
//  fragmentation CHIPS algorithm is the same in both models, but the first
//  interaction is temporary taken from the QGSC algorithm, which is not expected
//  to work at low energies, but temporary can be used as an interface to CHIPS.
//-------------------------------------------------------------------------------
#ifndef TQGSC_CHIPS_h
#define TQGSC_CHIPS_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TQGSC_CHIPS: public T
{
public:
  TQGSC_CHIPS(G4int ver = 1);
  virtual ~TQGSC_CHIPS();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};
#include "QGSC_CHIPS.icc"
typedef TQGSC_CHIPS<G4VModularPhysicsList> QGSC_CHIPS;

#endif
