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
// *                                                                  *
// * Parts of this code which have been  developed by Abdel-Waged     *
// * et al under contract (31-465) to the King Abdul-Aziz City for    *
// * Science and Technology (KACST), the National Centre of           *
// * Mathematics and Physics (NCMP), Saudi Arabia.                    *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/include/G4CRMCModel.hh
/// \brief Definition of the G4CRMCModel class
//
//---------------------------------------------------------------------------
//
// ClassName: G4CRMCModel
//
// Author:    2018 Alberto Ribon
//
// Geant4 wrapper hadronic model around CRMC.
// By changing a single number (model) it is possible to use one of the
// following 4 CRMC generators:  0 :  EPOS LHC
//                               1 :  EPOS 1.99
//                               6 :  SIBYLL 2.3c
//                              12 :  DPMJET 3
//
// We are pleased to acknowledge the contribution of Andrii Tykhonov
// (Universite' de Geneve, and member of the DAMPE Collaboration) 
// who made the first prototype interface between Geant4 and CRMC,
// and shared his expertise and implementation with us. 
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4CRMCModel_hh
#define G4CRMCModel_hh

#include "G4Nucleus.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadFinalState.hh"

#include <fstream>                
#include <string>

class G4IonTable;

#ifdef G4_USE_CRMC
  #include "CRMCinterface.h"
  extern CRMCdata gCRMC_data;
#else
  class CRMCinterface;
#endif


class G4CRMCModel : public G4HadronicInteraction {
  public:
    // ***LOOKHERE*** Select default model: 
    //                EPOS LHC (0), EPOS 1.99 (1), SIBYLL 2.3c (6), DPMJET 3 (12)
    G4CRMCModel( const G4int model = 0 );  
    virtual ~G4CRMCModel ();
    G4HadFinalState* ApplyYourself( const G4HadProjectile &theProjectile, G4Nucleus &theNucleus );
    G4ParticleDefinition* GetParticleDefinition( long particle_id );

    // The following virtual method can be used to override the default check levels
    // (relative energy violation value of 2% and absolute energy violation value of 1 GeV)
    // used to check (by the hadronic process, in the method G4HadronicProcess::CheckResult)
    // whether the (inelastic) final state is acceptable or not in terms of energy
    // conservation. If it is not, then the final-state is rejected, a "JustWarning"
    // exception is thrown, and another attempt is done to find an acceptable final-state.
    // Only in the case of 100 consecutive failed attempts, the program crashes.
    // In the case of EPOS-LHC, there are some cases of final-states that are rejected
    // because they do not pass the default energy conservation. The number of these
    // cases is relatively low and harmless (because these final-states are rejected),
    // therefore we suggest to keep using the default. If you want instead to keep
    // more final-states, i.e. be more tolerant on the energy violations, then uncomment
    // the following line and its implementation (in the source file).
    //virtual const std::pair< G4double, G4double > GetFatalEnergyCheckLevels() const;

  private: 
    G4bool operator==( G4CRMCModel& right );
    G4bool operator!=( G4CRMCModel& right );
    void WelcomeMessage () const;
    G4int CurrentEvent;
    G4int verbose;
    G4int fModel;
    G4HadFinalState  fFinalState;
    CRMCinterface*   fInterface;
    G4ParticleTable* fParticleTable;
    G4IonTable*      fIonTable;
};

#endif

