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
 *       Filename:  CexmcProductionModel.hh
 *
 *    Description:  interface to production model
 *
 *        Version:  1.0
 *        Created:  03.11.2009 16:50:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PRODUCTION_MODEL_HH
#define CEXMC_PRODUCTION_MODEL_HH

#include <G4Types.hh>
#include <G4String.hh>
#include <G4ios.hh>
#include <G4ParticleDefinition.hh>
#include "CexmcAngularRange.hh"
#include "CexmcProductionModelData.hh"
#include "CexmcHistoManager.hh"
#include "CexmcException.hh"
#include "CexmcCommon.hh"

class  CexmcProductionModelMessenger;


class  CexmcProductionModel
{
    public:
        explicit CexmcProductionModel( const G4String &  name = "unspecified",
                                       G4bool  fermiMotionIsOn = false );

        virtual ~CexmcProductionModel();

    public:
        void  ApplyFermiMotion( G4bool  on, G4bool  fromMessenger = true );

        void  SetAngularRange( G4double  top, G4double  bottom,
                               G4int  nmbOfDivs );

        void  SetAngularRanges( const CexmcAngularRangeList &  angularRanges_ );

        void  AddAngularRange( G4double  top, G4double  bottom,
                               G4int  nmbOfDivs );

        void  SetProductionModelData(
                    const CexmcProductionModelData &  productionModelData_ );

        void  PrintInitialData( void ) const;

        const CexmcAngularRangeList &  GetAngularRanges( void ) const;

        const CexmcAngularRangeList &  GetTriggeredAngularRanges( void ) const;

        const CexmcProductionModelData &  GetProductionModelData( void ) const;

        G4bool  IsFermiMotionOn( void ) const;

        void    SetTriggeredAngularRanges( G4double  opCosThetaSCM );

        const G4String &  GetName( void ) const;

    public:
        G4ParticleDefinition *  GetIncidentParticle( void ) const;

        G4ParticleDefinition *  GetNucleusParticle( void ) const;

        G4ParticleDefinition *  GetOutputParticle( void ) const;

        G4ParticleDefinition *  GetNucleusOutputParticle( void ) const;

    protected:
        virtual void            FermiMotionStatusChangeHook( void );

    private:
        G4bool  IsValidCandidateForAngularRange( G4double  top,
                             G4double  bottom, G4int  nmbOfDivs = 1 ) const;

        G4bool  IsGoodCandidateForAngularRange( G4double  top,
                                                G4double  bottom ) const;

    protected:
        G4String                  name;

        G4bool                    fermiMotionIsOn;

        CexmcAngularRangeList     angularRanges;

        CexmcAngularRangeList     angularRangesRef;

        CexmcAngularRangeList     triggeredAngularRanges;

        CexmcProductionModelData  productionModelData;

    protected:
        G4ParticleDefinition *    incidentParticle;

        G4ParticleDefinition *    nucleusParticle;

        G4ParticleDefinition *    outputParticle;

        G4ParticleDefinition *    nucleusOutputParticle;

    private:
        CexmcProductionModelMessenger *  messenger;
};


inline void  CexmcProductionModel::ApplyFermiMotion( G4bool  on,
                                                     G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    fermiMotionIsOn = on;

    FermiMotionStatusChangeHook();
}


inline void  CexmcProductionModel::SetAngularRanges(
                                const CexmcAngularRangeList &  angularRanges_ )
{
    angularRangesRef = angularRanges_;
    angularRanges = angularRangesRef;
#ifdef CEXMC_USE_ROOT
    CexmcHistoManager::Instance()->SetupARHistos( angularRanges );
#endif
}


inline void  CexmcProductionModel::SetProductionModelData(
                        const CexmcProductionModelData &  productionModelData_ )
{
    productionModelData = productionModelData_;
}


inline void  CexmcProductionModel::PrintInitialData( void ) const
{
    const char *  fermiMotionMsg( "Fermi motion in the target is off" );
    if ( fermiMotionIsOn )
        fermiMotionMsg = "Fermi motion in the target is on";

    G4cout << CEXMC_LINE_START << fermiMotionMsg << G4endl;
    G4cout << CEXMC_LINE_START << "Angular ranges:" << angularRanges;
}


inline const CexmcAngularRangeList &
                CexmcProductionModel::GetAngularRanges( void ) const
{
    return angularRanges;
}


inline const CexmcAngularRangeList &
                CexmcProductionModel::GetTriggeredAngularRanges( void ) const
{
    return triggeredAngularRanges;
}


inline const CexmcProductionModelData &
                CexmcProductionModel::GetProductionModelData( void ) const
{
    return productionModelData;
}


inline G4bool  CexmcProductionModel::IsFermiMotionOn( void ) const
{
    return fermiMotionIsOn;
}


inline const G4String &  CexmcProductionModel::GetName( void ) const
{
    return name;
}


inline  G4ParticleDefinition *  CexmcProductionModel::GetIncidentParticle(
                                                                    void ) const
{
    return incidentParticle;
}


inline  G4ParticleDefinition *  CexmcProductionModel::GetNucleusParticle( void )
                                                                        const
{
    return nucleusParticle;
}


inline  G4ParticleDefinition *  CexmcProductionModel::GetOutputParticle( void )
                                                                        const
{
    return outputParticle;
}


inline  G4ParticleDefinition *  CexmcProductionModel::GetNucleusOutputParticle(
                                                                    void ) const
{
    return nucleusOutputParticle;
}


inline G4bool  CexmcProductionModel::IsValidCandidateForAngularRange(
                    G4double  top, G4double  bottom, G4int  nmbOfDivs ) const
{
    return top > bottom && top <= 1.0 && top > -1.0 && bottom < 1.0 &&
           bottom >= -1.0 && nmbOfDivs >= 1;
}


#endif

