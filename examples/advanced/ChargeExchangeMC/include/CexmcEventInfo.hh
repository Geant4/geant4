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
 *       Filename:  CexmcEventInfo.hh
 *
 *    Description:  event information passed to run manager
 *
 *        Version:  1.0
 *        Created:  04.12.2009 14:47:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_EVENT_INFO_HH
#define CEXMC_EVENT_INFO_HH

#include <G4Types.hh>
#include "G4VUserEventInformation.hh"


class  CexmcEventInfo : public G4VUserEventInformation
{
    public:
        CexmcEventInfo( G4bool  edTriggerIsOk, G4bool  tpTriggerIsOk,
                        G4bool  reconstructionIsOk );

    public:
        void    Print( void ) const;

    public:
        G4bool  EdTriggerIsOk( void ) const;

        G4bool  TpTriggerIsOk( void ) const;

        G4bool  ReconstructionIsOk( void ) const;

    private:
        G4bool  edTriggerIsOk;

        G4bool  tpTriggerIsOk;

        G4bool  reconstructionIsOk;
};


inline G4bool  CexmcEventInfo::EdTriggerIsOk( void ) const
{
    return edTriggerIsOk;
}


inline G4bool  CexmcEventInfo::TpTriggerIsOk( void ) const
{
    return tpTriggerIsOk;
}


inline G4bool  CexmcEventInfo::ReconstructionIsOk( void ) const
{
    return reconstructionIsOk;
}


#endif

