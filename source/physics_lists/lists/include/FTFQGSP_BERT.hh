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
// $Id: FTFQGSP_BERT.hh 66241 2012-12-13 18:34:42Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:  FTFQGSP_BERT
//
// Author: Alberto Ribon, October 2017
//   created from FTFP_BERT, using QGS string fragmentation instead of the
//   Lund string fragmentation
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef TFTFQGSP_BERT_h
#define TFTFQGSP_BERT_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VModularPhysicsList.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TFTFQGSP_BERT: public T
{
public:
  TFTFQGSP_BERT(G4int ver = 1);
  virtual ~TFTFQGSP_BERT();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};
#include "FTFQGSP_BERT.icc"
typedef TFTFQGSP_BERT<G4VModularPhysicsList> FTFQGSP_BERT;

#endif



