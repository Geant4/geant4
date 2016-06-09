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
/*
 * =============================================================================
 *
 *       Filename:  CexmcHadronicProcess.hh
 *
 *    Description:  hadronic process with production model
 *
 *        Version:  1.0
 *        Created:  31.10.2009 23:44:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_HADRONIC_PROCESS_HH
#define CEXMC_HADRONIC_PROCESS_HH

#include <G4HadronicProcess.hh>
#include <G4Nucleus.hh>
#include "CexmcCommon.hh"

class  G4VParticleChange;
class  G4ParticleDefinition;
class  G4Track;
class  G4Step;
class  G4Material;
class  G4HadronicInteraction;
class  CexmcProductionModel;


class  CexmcHadronicProcess : public G4HadronicProcess
{
    public:
        explicit CexmcHadronicProcess(
                        const G4String &  name = CexmcStudiedProcessLastName );

        ~CexmcHadronicProcess();

    public:
        G4VParticleChange *  PostStepDoIt( const G4Track &  track,
                                           const G4Step &  step );

        G4bool  IsApplicable( const G4ParticleDefinition &  particle );

    public:
        void    RegisterProductionModel( CexmcProductionModel *  model );

    private:
        void    CalculateTargetNucleus( const G4Material *  material );

        void    FillTotalResult( G4HadFinalState *  hadFinalState,
                                 const G4Track &  track );

    private:
        CexmcProductionModel *   productionModel;

        G4HadronicInteraction *  interaction;

    private:
        G4ParticleChange *      theTotalResult;

        G4Nucleus               targetNucleus;

        G4bool                  isInitialized;
};


#endif

