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
/// \file hadronic/Hadr02/include/HIJING.hh
/// \brief Definition of the HIJING class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2012  Andrea Dotti
//   created from FTFP_BERT
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef THIJING_h
#define THIJING_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class THIJING: public T
{
public:
  THIJING(G4int ver = 1);
  virtual ~THIJING();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};
#ifdef G4_USE_HIJING
#include "HIJING.icc"
#else
template<class T>
THIJING<T>::THIJING(G4int) : T()
{
  G4ExceptionDescription de;
  de<<"Support for HIJING not enabled"<<G4endl;
  G4Exception(__FILE__,"HIJING-01",FatalException,de,
  "Code should be compiled with G4_USE_HIJING environment variable set.");
}

template<class T>
THIJING<T>::~THIJING() { }
template<class T>
void THIJING<T>::SetCuts() { }
#endif

typedef THIJING<G4VModularPhysicsList> HIJING;

#endif //THIJING_h


