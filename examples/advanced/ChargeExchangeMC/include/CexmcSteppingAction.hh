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
 * ============================================================================
 *
 *       Filename:  CexmcSteppingAction.hh
 *
 *    Description:  stepping action
 *
 *        Version:  1.0
 *        Created:  27.10.2009 15:59:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#ifndef CEXMC_STEPPING_ACTION_HH
#define CEXMC_STEPPING_ACTION_HH

#include <G4UserSteppingAction.hh>

class  G4Step;
class  G4LogicalVolume;
class  CexmcPhysicsManager;


class  CexmcSteppingAction : public G4UserSteppingAction
{
    public:
        explicit CexmcSteppingAction( CexmcPhysicsManager *  physicsManager );

    public:
        void  UserSteppingAction( const G4Step *  step );

    private:
        CexmcPhysicsManager *    physicsManager;

        const G4LogicalVolume *  targetVolume;
};


#endif

