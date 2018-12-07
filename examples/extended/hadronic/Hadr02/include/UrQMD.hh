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
/// \file hadronic/Hadr02/include/UrQMD.hh
/// \brief Definition of the UrQMD class
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
#ifndef TUrQMD_h
#define TUrQMD_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TUrQMD: public T
{
public:
  TUrQMD(G4int ver = 1);
  virtual ~TUrQMD();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};
#ifdef G4_USE_URQMD
#include "UrQMD.icc"
#else
template<class T>
TUrQMD<T>::TUrQMD(G4int) : T()
{
  G4ExceptionDescription de;
  de<<"Support for UrQMD not enabled"<<G4endl;
  G4Exception(__FILE__,"UrQMD-01",FatalException,de,
    "Code should be compiled with G4_USE_URQMD environment variable set.");
}

template<class T>
TUrQMD<T>::~TUrQMD() { }
template<class T>
void TUrQMD<T>::SetCuts() { }
#endif

typedef TUrQMD<G4VModularPhysicsList> UrQMD;

#endif //TUrQMD_h


