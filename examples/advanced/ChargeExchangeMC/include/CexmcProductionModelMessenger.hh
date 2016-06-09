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
 *       Filename:  CexmcProductionModelMessenger.hh
 *
 *    Description:  set various production model aspects
 *
 *        Version:  1.0
 *        Created:  03.11.2009 15:50:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PRODUCTION_MODEL_MESSENGER_HH
#define CEXMC_PRODUCTION_MODEL_MESSENGER_HH

#include <G4UImessenger.hh>

class  G4UIcommand;
class  G4UIcmdWithABool;
class  G4UIcmdWith3Vector;
class  CexmcProductionModel;


class  CexmcProductionModelMessenger : public G4UImessenger
{
    public:
        explicit CexmcProductionModelMessenger( CexmcProductionModel *
                                                productionModel );

        ~CexmcProductionModelMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcProductionModel *  productionModel;

        G4UIcmdWithABool *      applyFermiMotion;

        G4UIcmdWith3Vector *    setAngularRange;

        G4UIcmdWith3Vector *    addAngularRange;
};


#endif

