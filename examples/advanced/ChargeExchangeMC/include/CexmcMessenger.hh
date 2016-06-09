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
 *       Filename:  CexmcMessenger.hh
 *
 *    Description:  common messenger stuff (directories etc.)
 *
 *        Version:  1.0
 *        Created:  15.11.2009 12:48:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_MESSENGER_HH
#define CEXMC_MESSENGER_HH

#include <G4String.hh>

class  G4UIdirectory;


class  CexmcMessenger
{
    public:
        static CexmcMessenger *  Instance( void );

        static void              Destroy( void );

    private:
        CexmcMessenger();

        ~CexmcMessenger();

    public:
        static G4String  mainDirName;

        static G4String  geometryDirName;

        static G4String  physicsDirName;

        static G4String  gunDirName;

        static G4String  detectorDirName;

        static G4String  eventDirName;

        static G4String  runDirName;

        static G4String  monitorDirName;

        static G4String  targetDirName;

        static G4String  vetoCounterDirName;

        static G4String  vetoCounterLeftDirName;

        static G4String  vetoCounterRightDirName;

        static G4String  calorimeterDirName;

        static G4String  calorimeterLeftDirName;

        static G4String  calorimeterRightDirName;

        static G4String  monitorEDDirName;

        static G4String  vetoCounterEDDirName;

        static G4String  vetoCounterLeftEDDirName;

        static G4String  vetoCounterRightEDDirName;

        static G4String  calorimeterEDDirName;

        static G4String  calorimeterLeftEDDirName;

        static G4String  calorimeterRightEDDirName;

        static G4String  reconstructorDirName;

        static G4String  visDirName;

#ifdef CEXMC_USE_ROOT
        static G4String  histoDirName;
#endif

    private:
        static CexmcMessenger *  instance;

    private:
        G4UIdirectory *  mainDir;

        G4UIdirectory *  geometryDir;

        G4UIdirectory *  physicsDir;

        G4UIdirectory *  gunDir;

        G4UIdirectory *  detectorDir;

        G4UIdirectory *  eventDir;

        G4UIdirectory *  runDir;

        G4UIdirectory *  monitorDir;

        G4UIdirectory *  targetDir;

        G4UIdirectory *  vetoCounterDir;

        G4UIdirectory *  vetoCounterLeftDir;

        G4UIdirectory *  vetoCounterRightDir;

        G4UIdirectory *  calorimeterDir;

        G4UIdirectory *  calorimeterLeftDir;

        G4UIdirectory *  calorimeterRightDir;

        G4UIdirectory *  monitorEDDir;

        G4UIdirectory *  vetoCounterEDDir;

        G4UIdirectory *  vetoCounterLeftEDDir;

        G4UIdirectory *  vetoCounterRightEDDir;

        G4UIdirectory *  calorimeterEDDir;

        G4UIdirectory *  calorimeterLeftEDDir;

        G4UIdirectory *  calorimeterRightEDDir;

        G4UIdirectory *  reconstructorDir;

        G4UIdirectory *  visDir;

#ifdef CEXMC_USE_ROOT
        G4UIdirectory *  histoDir;
#endif
};


#endif

