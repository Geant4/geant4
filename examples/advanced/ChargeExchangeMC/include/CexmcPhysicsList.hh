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
 *       Filename:  CexmcPhysicsList.hh
 *
 *    Description:  mandatory physics list
 *
 *        Version:  1.0
 *        Created:  11.10.2009 14:51:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PHYSICS_LIST_HH
#define CEXMC_PHYSICS_LIST_HH

#include <Randomize.hh>
#include <G4Track.hh>
#include <G4StepPoint.hh>
#include <G4ThreeVector.hh>
#include <G4AffineTransform.hh>
#include "CexmcPhysicsManager.hh"
#include "CexmcProductionModel.hh"
#include "CexmcIncidentParticleTrackInfo.hh"
#include "CexmcSetup.hh"
#include "CexmcException.hh"
#include "CexmcCommon.hh"


template  < typename  BasePhysics, template  < typename > class  StudiedPhysics,
            typename  ProductionModel >
class  CexmcPhysicsList : public BasePhysics, public CexmcPhysicsManager
{
    public:
        CexmcPhysicsList();

    public:
        CexmcProductionModel *  GetProductionModel( void );

        G4bool  IsStudiedProcessAllowed( void ) const;

        void    ResampleTrackLengthInTarget( const G4Track *  track,
                                             const G4StepPoint *  stepPoint );

        void    SetupConstructionHook( const CexmcSetup *  setup );

    protected:
        void  CalculateBasicMaxIL( const G4ThreeVector &  direction );

    private:
        StudiedPhysics< ProductionModel > *  studiedPhysics;

        G4VSolid *                           targetSolid;

        G4AffineTransform                    targetTransform;
};


template  < typename  BasePhysics, template  < typename > class  StudiedPhysics,
            typename  ProductionModel >
CexmcPhysicsList< BasePhysics, StudiedPhysics, ProductionModel >::
                CexmcPhysicsList() : studiedPhysics( NULL ), targetSolid( NULL )
{
    studiedPhysics = new StudiedPhysics< ProductionModel >( this );
    this->RegisterPhysics( studiedPhysics );
}


template  < typename  BasePhysics, template  < typename > class  StudiedPhysics,
            typename  ProductionModel >
CexmcProductionModel *
            CexmcPhysicsList< BasePhysics, StudiedPhysics, ProductionModel >::
                GetProductionModel( void )
{
    return studiedPhysics->GetProductionModel();
}


template  < typename  BasePhysics, template  < typename > class  StudiedPhysics,
            typename  ProductionModel >
G4bool  CexmcPhysicsList< BasePhysics, StudiedPhysics, ProductionModel >::
                IsStudiedProcessAllowed( void ) const
{
    return numberOfTriggeredStudiedInteractions == 0;
}


template  < typename  BasePhysics, template  < typename > class  StudiedPhysics,
            typename  ProductionModel >
void  CexmcPhysicsList< BasePhysics, StudiedPhysics, ProductionModel >::
                ResampleTrackLengthInTarget( const G4Track *  track,
                                             const G4StepPoint *  stepPoint )
{
    /* BEWARE: all callers must ensure that:
     * 1) track (or stepPoint if not NULL) is inside target volume:
     *    in this case we can use already calculated targetTransform
     * 2) track info object is of type CexmcIncidentParticleTrackInfo*:
     *    in this case we can use static_cast<> for trackInfo */
    CexmcIncidentParticleTrackInfo *  trackInfo(
                static_cast< CexmcIncidentParticleTrackInfo * >(
                                                track->GetUserInformation() ) );

    if ( ! trackInfo )
        return;

    G4ThreeVector  position;
    G4ThreeVector  direction;

    if ( stepPoint )
    {
        position = targetTransform.TransformPoint( stepPoint->GetPosition() );
        direction = targetTransform.TransformAxis(
                                            stepPoint->GetMomentumDirection() );
    }
    else
    {
        position = targetTransform.TransformPoint( track->GetPosition() );
        direction = targetTransform.TransformAxis(
                                            track->GetMomentumDirection() );
    }

    G4double  distanceInTarget( targetSolid->DistanceToOut( position,
                                                            direction ) );
    trackInfo->ResetCurrentTrackLengthInTarget();
    trackInfo->SetFinalTrackLengthInTarget( G4UniformRand() *
                                std::max( distanceInTarget, proposedMaxIL ) );
    trackInfo->SetNeedsTrackLengthResampling( false );
}


template  < typename  BasePhysics, template  < typename > class  StudiedPhysics,
            typename  ProductionModel >
void  CexmcPhysicsList< BasePhysics, StudiedPhysics, ProductionModel >::
                CalculateBasicMaxIL( const G4ThreeVector &  direction )
{
    /* basicMaxIL is double distance from the point (0, 0, 0) to the edge of the
     * target solid along the specified direction */
    basicMaxIL = targetSolid->DistanceToOut( G4ThreeVector(),
                            targetTransform.TransformAxis( direction ) ) * 2;
}


template  < typename  BasePhysics, template  < typename > class  StudiedPhysics,
            typename  ProductionModel >
void  CexmcPhysicsList< BasePhysics, StudiedPhysics, ProductionModel >::
                SetupConstructionHook( const CexmcSetup *  setup )
{
    targetSolid = setup->GetVolume( CexmcSetup::Target )->GetSolid();
    targetTransform = setup->GetTargetTransform().Inverse();
}


#endif

