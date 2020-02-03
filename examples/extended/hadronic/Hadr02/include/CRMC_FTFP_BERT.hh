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
/// \file hadronic/Hadr02/include/CRMC_FTFP_BERT.hh
/// \brief Definition of the CRMC_FTFP_BERT class
//
//
//---------------------------------------------------------------------------
//
// ClassName: CRMC_FTFP_BERT  
//
// Author:    2018 Alberto Ribon
//
// This is a variant of the FTFP_BERT physics list, whereby CRMC is used 
// for modeling final-state pion- , kaon- , proton- , neutron- and 
// ion-nuclear inelastic interactions at very high energies (above a
// threshold defined in HadronPhysicsCRMC_FTFP_BERT and IonCRMCPhysics).
// For the remaining inelastic interactions (i.e. hyperon- , antihyperon- ,
// antinucleon- and light anti-ion-nuclear interactions, as well as
// for all elastic final-state interactions, and for all elastic and
// inelastic hadronic cross sections, the usual Geant4 approach, as in
// FTFP_BERT is used.
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef TCRMC_FTFP_BERT_h
#define TCRMC_FTFP_BERT_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"


template<class T> 
class TCRMC_FTFP_BERT: public T {
  public:
    TCRMC_FTFP_BERT( G4int ver = 1 );
    virtual ~TCRMC_FTFP_BERT();
    virtual void SetCuts();
  private:
    enum { ok = CompileTimeConstraints::IsA< T, G4VModularPhysicsList >::ok };
};


#ifdef G4_USE_CRMC
#include "CRMC_FTFP_BERT.icc"
#else
template<class T>
TCRMC_FTFP_BERT<T>::TCRMC_FTFP_BERT( G4int ) : T() {
  G4ExceptionDescription de;
  de << "Support for CRMC_FTFP_BERT not enabled" << G4endl;
  G4Exception( __FILE__, "CRMC_FTFP_BERT-01", FatalException, de,
              "Code should be compiled with G4_USE_CRMC environment variable set.");
}

template<class T>
TCRMC_FTFP_BERT<T>::~TCRMC_FTFP_BERT() { }
template<class T>
void TCRMC_FTFP_BERT<T>::SetCuts() {}
#endif


typedef TCRMC_FTFP_BERT< G4VModularPhysicsList > CRMC_FTFP_BERT;

#endif //TCRMC_FTFP_BERT_h

