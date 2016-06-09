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
 *       Filename:  CexmcStudiedProcess.hh
 *
 *    Description:  studied process in the target
 *
 *        Version:  1.0
 *        Created:  26.10.2009 20:41:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_STUDIED_PROCESS_HH
#define CEXMC_STUDIED_PROCESS_HH

#include <G4WrapperProcess.hh>
#include <G4ProcessType.hh>

class  G4VParticleChange;
class  CexmcPhysicsManager;


class  CexmcStudiedProcess : public G4WrapperProcess
{
    public:
        explicit  CexmcStudiedProcess( CexmcPhysicsManager *  physicsManager,
                                   G4ProcessType  processType = fUserDefined );

    public:
        G4double  PostStepGetPhysicalInteractionLength( const G4Track &  track,
                    G4double  previousStepSize, G4ForceCondition *  condition );

        G4VParticleChange *  PostStepDoIt( const G4Track &  track,
                                           const G4Step &  step );

    private:
        CexmcPhysicsManager *  physicsManager;
};


#endif

