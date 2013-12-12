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
 *       Filename:  CexmcEnergyDepositDigitizer.hh
 *
 *    Description:  digitizes of energy deposit in a single event
 *
 *        Version:  1.0
 *        Created:  23.11.2009 14:14:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_ENERGY_DEPOSIT_DIGITIZER_HH
#define CEXMC_ENERGY_DEPOSIT_DIGITIZER_HH

#include <iosfwd>
#include <G4VDigitizerModule.hh>
#include <G4SystemOfUnits.hh>
#include "CexmcEnergyDepositStore.hh"
#include "CexmcSimpleRangeWithValue.hh"
#include "CexmcException.hh"
#include "CexmcCommon.hh"

class  G4String;
class  CexmcEnergyDepositDigitizerMessenger;


class  CexmcEnergyDepositDigitizer : public G4VDigitizerModule
{
    public:
        explicit CexmcEnergyDepositDigitizer( const G4String &  name );

        ~CexmcEnergyDepositDigitizer();

    public:
        void      Digitize( void );

    public:
        G4double  GetMonitorED( void ) const;

        G4double  GetVetoCounterEDLeft( void ) const;

        G4double  GetVetoCounterEDRight( void ) const;

        G4double  GetCalorimeterEDLeft( void ) const;

        G4double  GetCalorimeterEDRight( void ) const;

        G4int     GetCalorimeterEDLeftMaxX( void ) const;

        G4int     GetCalorimeterEDLeftMaxY( void ) const;

        G4int     GetCalorimeterEDRightMaxX( void ) const;

        G4int     GetCalorimeterEDRightMaxY( void ) const;

        const CexmcEnergyDepositCalorimeterCollection &
                                GetCalorimeterEDLeftCollection( void ) const;

        const CexmcEnergyDepositCalorimeterCollection &
                                GetCalorimeterEDRightCollection( void ) const;

    public:
        G4bool    MonitorHasTriggered( void ) const;

        G4bool    HasTriggered( void ) const;

    public:
        void      SetMonitorThreshold( G4double  value,
                                       G4bool  fromMessenger = true );

        void      SetVetoCounterLeftThreshold( G4double  value,
                                               G4bool  fromMessenger = true );

        void      SetVetoCounterRightThreshold( G4double  value,
                                                G4bool  fromMessenger = true );

        void      SetVetoCountersThreshold( G4double  value );

        void      SetCalorimeterLeftThreshold( G4double  value,
                                               G4bool  fromMessenger = true );

        void      SetCalorimeterRightThreshold( G4double  value,
                                                G4bool  fromMessenger = true );

        void      SetCalorimetersThreshold( G4double  value );

        void      SetCalorimeterTriggerAlgorithm(
                                    CexmcCalorimeterTriggerAlgorithm  value,
                                    G4bool  fromMessenger = true );

        void      SetOuterCrystalsVetoAlgorithm(
                                    CexmcOuterCrystalsVetoAlgorithm  value,
                                    G4bool  fromMessenger = true );

        void      SetOuterCrystalsVetoFraction( G4double  value,
                                                G4bool  fromMessenger = true );

        void      ApplyFiniteCrystalResolution( G4bool  value,
                                                G4bool  fromMessenger = true );

        void      AddCrystalResolutionRange( G4double  bottom, G4double  top,
                                             G4double  value,
                                             G4bool  fromMessenger = true );

        void      ClearCrystalResolutionData( G4bool  fromMessenger = true );

        void      SetCrystalResolutionData(
                            const CexmcEnergyRangeWithDoubleValueList &  data );

        G4double  GetMonitorThreshold( void ) const;

        G4double  GetVetoCounterLeftThreshold( void ) const;

        G4double  GetVetoCounterRightThreshold( void ) const;

        G4double  GetCalorimeterLeftThreshold( void ) const;

        G4double  GetCalorimeterRightThreshold( void ) const;

        CexmcCalorimeterTriggerAlgorithm
                  GetCalorimeterTriggerAlgorithm( void ) const;

        CexmcOuterCrystalsVetoAlgorithm
                  GetOuterCrystalsVetoAlgorithm( void ) const;

        G4double  GetOuterCrystalsVetoFraction( void ) const;

        G4bool    IsFiniteCrystalResolutionApplied( void ) const;

        const CexmcEnergyRangeWithDoubleValueList &
                  GetCrystalResolutionData( void ) const;

    public:
        G4bool    IsOuterCrystal( G4int  column, G4int  row ) const;

    private:
        void      InitializeData( void );

    private:
        G4double                                 monitorED;

        G4double                                 vetoCounterEDLeft;

        G4double                                 vetoCounterEDRight;

        CexmcEnergyDepositCalorimeterCollection  calorimeterEDLeftCollection;

        CexmcEnergyDepositCalorimeterCollection  calorimeterEDRightCollection;

        G4double                                 calorimeterEDLeft;

        G4double                                 calorimeterEDRight;

        G4int                                    calorimeterEDLeftMaxX;

        G4int                                    calorimeterEDLeftMaxY;

        G4int                                    calorimeterEDRightMaxX;

        G4int                                    calorimeterEDRightMaxY;

        G4bool                                   monitorHasTriggered;

        G4bool                                   hasTriggered;

    private:
        G4double                                 monitorEDThreshold;

        G4double                                 vetoCounterEDLeftThreshold;

        G4double                                 vetoCounterEDRightThreshold;

        G4double                                 calorimeterEDLeftThreshold;

        G4double                                 calorimeterEDRightThreshold;

        CexmcCalorimeterTriggerAlgorithm         calorimeterTriggerAlgorithm;

        CexmcOuterCrystalsVetoAlgorithm          outerCrystalsVetoAlgorithm;

        G4double                                 outerCrystalsVetoFraction;

        G4double                                 monitorEDThresholdRef;

        G4double                                 vetoCounterEDLeftThresholdRef;

        G4double                                 vetoCounterEDRightThresholdRef;

        G4double                                 calorimeterEDLeftThresholdRef;

        G4double                                 calorimeterEDRightThresholdRef;

        CexmcCalorimeterTriggerAlgorithm         calorimeterTriggerAlgorithmRef;

        CexmcOuterCrystalsVetoAlgorithm          outerCrystalsVetoAlgorithmRef;

        G4double                                 outerCrystalsVetoFractionRef;

    private:
        G4int                                    nCrystalsInColumn;

        G4int                                    nCrystalsInRow;

    private:
        G4bool                                   applyFiniteCrystalResolution;

        CexmcEnergyRangeWithDoubleValueList      crystalResolutionData;

    private:
        CexmcEnergyDepositDigitizerMessenger *   messenger;
};


inline G4double  CexmcEnergyDepositDigitizer::GetMonitorED( void ) const
{
    return monitorED;
}


inline G4double  CexmcEnergyDepositDigitizer::GetVetoCounterEDLeft( void ) const
{
    return vetoCounterEDLeft;
}


inline G4double  CexmcEnergyDepositDigitizer::GetVetoCounterEDRight( void )
                                                                        const
{
    return vetoCounterEDRight;
}


inline G4double  CexmcEnergyDepositDigitizer::GetCalorimeterEDLeft( void ) const
{
    return calorimeterEDLeft;
}


inline G4double  CexmcEnergyDepositDigitizer::GetCalorimeterEDRight( void )
                                                                        const
{
    return calorimeterEDRight;
}


inline G4int  CexmcEnergyDepositDigitizer::GetCalorimeterEDLeftMaxX( void )
                                                                        const
{
    return calorimeterEDLeftMaxX;
}


inline G4int  CexmcEnergyDepositDigitizer::GetCalorimeterEDLeftMaxY( void )
                                                                        const
{
    return calorimeterEDLeftMaxY;
}


inline G4int  CexmcEnergyDepositDigitizer::GetCalorimeterEDRightMaxX( void )
                                                                        const
{
    return calorimeterEDRightMaxX;
}


inline G4int  CexmcEnergyDepositDigitizer::GetCalorimeterEDRightMaxY( void )
                                                                        const
{
    return calorimeterEDRightMaxY;
}


inline const CexmcEnergyDepositCalorimeterCollection &
    CexmcEnergyDepositDigitizer::GetCalorimeterEDLeftCollection( void ) const
{
    return calorimeterEDLeftCollection;
}


inline const CexmcEnergyDepositCalorimeterCollection &
    CexmcEnergyDepositDigitizer::GetCalorimeterEDRightCollection( void ) const
{
    return calorimeterEDRightCollection;
}


inline G4bool  CexmcEnergyDepositDigitizer::MonitorHasTriggered( void ) const
{
    return monitorHasTriggered;
}


inline G4bool  CexmcEnergyDepositDigitizer::HasTriggered( void ) const
{
    return hasTriggered;
}


inline void  CexmcEnergyDepositDigitizer::SetMonitorThreshold(
                                      G4double  value, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcBadThreshold,
                                       value < monitorEDThresholdRef );
    else
        monitorEDThresholdRef = value;

    monitorEDThreshold = value;
}


inline void  CexmcEnergyDepositDigitizer::SetVetoCounterLeftThreshold(
                                      G4double  value, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcBadThreshold,
                                       value > vetoCounterEDLeftThresholdRef );
    else
        vetoCounterEDLeftThresholdRef = value;

    vetoCounterEDLeftThreshold = value;
}


inline void  CexmcEnergyDepositDigitizer::SetVetoCounterRightThreshold(
                                      G4double  value, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcBadThreshold,
                                       value > vetoCounterEDRightThresholdRef );
    else
        vetoCounterEDRightThresholdRef = value;

    vetoCounterEDRightThreshold = value;
}


inline void  CexmcEnergyDepositDigitizer::SetVetoCountersThreshold(
                                                            G4double  value )
{
    ThrowExceptionIfProjectIsRead( CexmcBadThreshold,
                                   value > vetoCounterEDLeftThresholdRef ||
                                   value > vetoCounterEDRightThresholdRef );

    vetoCounterEDLeftThreshold = value;
    vetoCounterEDRightThreshold = value;
}


inline void  CexmcEnergyDepositDigitizer::SetCalorimeterLeftThreshold(
                                      G4double  value, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcBadThreshold,
                                       value < calorimeterEDLeftThresholdRef );
    else
        calorimeterEDLeftThresholdRef = value;

    calorimeterEDLeftThreshold = value;
}


inline void  CexmcEnergyDepositDigitizer::SetCalorimeterRightThreshold(
                                      G4double  value, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcBadThreshold,
                                       value < calorimeterEDRightThresholdRef );
    else
        calorimeterEDRightThresholdRef = value;

    calorimeterEDRightThreshold = value;
}


inline void  CexmcEnergyDepositDigitizer::SetCalorimetersThreshold(
                                                            G4double  value )
{
    ThrowExceptionIfProjectIsRead( CexmcBadThreshold,
                                   value < calorimeterEDLeftThresholdRef ||
                                   value < calorimeterEDRightThresholdRef );

    calorimeterEDLeftThreshold = value;
    calorimeterEDRightThreshold = value;
}


inline void  CexmcEnergyDepositDigitizer::SetCalorimeterTriggerAlgorithm(
                CexmcCalorimeterTriggerAlgorithm  value, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcBadCalorimeterTriggerAlgorithm,
            ! ( calorimeterTriggerAlgorithmRef ==
                                       CexmcAllCrystalsMakeEDTriggerThreshold ||
                value == calorimeterTriggerAlgorithmRef ) );
    else
        calorimeterTriggerAlgorithmRef = value;

    calorimeterTriggerAlgorithm = value;
}


inline void  CexmcEnergyDepositDigitizer::SetOuterCrystalsVetoAlgorithm(
                CexmcOuterCrystalsVetoAlgorithm  value, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcBadOCVetoAlgorithm,
            ! ( outerCrystalsVetoAlgorithmRef == CexmcNoOuterCrystalsVeto ||
                value == outerCrystalsVetoAlgorithmRef ) );
    else
        outerCrystalsVetoAlgorithmRef = value;

    outerCrystalsVetoAlgorithm = value;
}


inline void  CexmcEnergyDepositDigitizer::SetOuterCrystalsVetoFraction(
                                      G4double  value, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcBadOCVetoFraction,
                                       value > outerCrystalsVetoFractionRef );
    else
        outerCrystalsVetoFractionRef = value;

    outerCrystalsVetoFraction = value;
}


inline void  CexmcEnergyDepositDigitizer::ApplyFiniteCrystalResolution(
                                      G4bool  value, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    applyFiniteCrystalResolution = value;
}


inline void  CexmcEnergyDepositDigitizer::AddCrystalResolutionRange(
                                      G4double  bottom, G4double  top,
                                      G4double  value, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    /* range boundaries are given in GeV */
    crystalResolutionData.push_back( CexmcEnergyRangeWithDoubleValue(
                                            bottom * GeV, top * GeV, value ) );
}


inline void  CexmcEnergyDepositDigitizer::ClearCrystalResolutionData(
                                                G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    crystalResolutionData.clear();
}


inline void  CexmcEnergyDepositDigitizer::SetCrystalResolutionData(
                            const CexmcEnergyRangeWithDoubleValueList &  data )
{
    ClearCrystalResolutionData( false );
    crystalResolutionData = data;
}


inline G4bool  CexmcEnergyDepositDigitizer::IsOuterCrystal( G4int  column,
                                                            G4int  row ) const
{
    return column == 0 || column == nCrystalsInRow - 1 ||
           row == 0 || row == nCrystalsInColumn - 1;
}


inline G4double  CexmcEnergyDepositDigitizer::GetMonitorThreshold( void ) const
{
    return monitorEDThreshold;
}


inline G4double  CexmcEnergyDepositDigitizer::GetVetoCounterLeftThreshold(
                                                                    void ) const
{
    return vetoCounterEDLeftThreshold;
}


inline G4double  CexmcEnergyDepositDigitizer::GetVetoCounterRightThreshold(
                                                                    void ) const
{
    return vetoCounterEDRightThreshold;
}


inline G4double  CexmcEnergyDepositDigitizer::GetCalorimeterLeftThreshold(
                                                                    void ) const
{
    return calorimeterEDLeftThreshold;
}


inline G4double  CexmcEnergyDepositDigitizer::GetCalorimeterRightThreshold(
                                                                    void ) const
{
    return calorimeterEDRightThreshold;
}


inline  CexmcCalorimeterTriggerAlgorithm
                CexmcEnergyDepositDigitizer::GetCalorimeterTriggerAlgorithm(
                                                                    void ) const
{
    return calorimeterTriggerAlgorithm;
}


inline  CexmcOuterCrystalsVetoAlgorithm
                CexmcEnergyDepositDigitizer::GetOuterCrystalsVetoAlgorithm(
                                                                    void ) const
{
    return outerCrystalsVetoAlgorithm;
}


inline G4double  CexmcEnergyDepositDigitizer::GetOuterCrystalsVetoFraction(
                                                                    void ) const
{
    return outerCrystalsVetoFraction;
}


inline G4bool  CexmcEnergyDepositDigitizer::IsFiniteCrystalResolutionApplied(
                                                                    void ) const
{
    return applyFiniteCrystalResolution;
}


inline const CexmcEnergyRangeWithDoubleValueList &
            CexmcEnergyDepositDigitizer::GetCrystalResolutionData( void ) const
{
    return crystalResolutionData;
}


std::ostream &  operator<<( std::ostream &  out,
                const CexmcEnergyDepositCalorimeterCollection &  edCollection );


#endif

