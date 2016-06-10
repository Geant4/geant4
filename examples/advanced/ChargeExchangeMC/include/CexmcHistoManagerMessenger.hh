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
 *       Filename:  CexmcHistoManagerMessenger.hh
 *
 *    Description:  commands to list and show histograms
 *
 *        Version:  1.0
 *        Created:  17.12.2009 21:38:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_HISTO_MANAGER_MESSENGER_HH
#define CEXMC_HISTO_MANAGER_MESSENGER_HH

#ifdef CEXMC_USE_ROOT

#include <G4UImessenger.hh>

class  G4UIcommand;
class  G4UIcmdWithoutParameter;
class  G4UIcmdWithAnInteger;
class  G4UIcmdWithAString;
class  CexmcHistoManager;


class  CexmcHistoManagerMessenger : public G4UImessenger
{
    public:
        explicit CexmcHistoManagerMessenger(
                                        CexmcHistoManager *  histoManager );

        ~CexmcHistoManagerMessenger();

    public:
        void  SetNewValue( G4UIcommand *  cmd, G4String  value );

    private:
        CexmcHistoManager *        histoManager;

        G4UIcmdWithAnInteger *     setVerboseLevel;

        G4UIcmdWithoutParameter *  listHistos;

        G4UIcmdWithAString *       printHisto;

#ifdef CEXMC_USE_ROOTQT
        G4UIcmdWithAString *       drawHisto;

        G4UIcmdWithAString *       addHistoMenu;

        G4UIcmdWithAString *       drawHistoOptions1D;

        G4UIcmdWithAString *       drawHistoOptions2D;

        G4UIcmdWithAString *       drawHistoOptions3D;
#endif
};

#endif

#endif

