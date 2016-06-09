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
 *       Filename:  CexmcScenePrimitives.hh
 *
 *    Description:  auxiliary scene primitives (radial lines etc.)
 *
 *        Version:  1.0
 *        Created:  03.01.2011 11:27:34
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SCENE_PRIMITIVES_HH
#define CEXMC_SCENE_PRIMITIVES_HH

#include <vector>
#include <map>
#include <G4Colour.hh>
#include <G4ThreeVector.hh>
#include <G4VModel.hh>
#include <G4VVisManager.hh>

class  G4VGraphicsScene;
class  CexmcSetup;
class  CexmcScenePrimitivesMessenger;


enum  CexmcSPType
{
    CexmcTargetCenterMark_SP,
    CexmcRadialLine_SP,
    CexmcInnerCrystalsHl_SP
};


class  CexmcScenePrimitives : public G4VModel
{
    private:
        struct  CexmcRadialLine
        {
            CexmcRadialLine( const G4ThreeVector &  line ) :
                theta( line.x() ), phi( line.y() ), length( line.z() )
            {}

            G4double  theta;

            G4double  phi;

            G4double  length;
        };

        typedef std::vector< CexmcRadialLine >      CexmcRadialLines;

        typedef std::map< CexmcSPType, G4Colour >   CexmcSPColourMap;

    public:
        explicit CexmcScenePrimitives( CexmcSetup *  setup );

        ~CexmcScenePrimitives();

    public:
        void  DescribeYourselfTo( G4VGraphicsScene &  scene );

    public:
        void  MarkTargetCenter( G4bool  on = true );

        void  DrawRadialLine( const G4ThreeVector &  line );

        void  HighlightInnerCrystals( G4bool = true );

        void  ClearRadialLines( void );

        void  SetColour( CexmcSPType  primitive, const G4Colour &  colour );

    private:
        void  DrawRadialLine( G4VGraphicsScene &  scene,
                              const CexmcRadialLine *  rLine );

        void  MarkTargetCenter( G4VGraphicsScene &  scene );

        void  HighlightInnerCrystals( G4VGraphicsScene &  scene );

    private:
        void  UpdateScene( void );

    private:
        CexmcSetup *                     setup;

        G4bool                           markTargetCenter;

        G4bool                           highlightInnerCrystals;

        CexmcRadialLines                 radialLines;

        CexmcSPColourMap                 spColours;

    private:
        CexmcScenePrimitivesMessenger *  messenger;
};


inline void  CexmcScenePrimitives::SetColour( CexmcSPType  primitive,
                                              const G4Colour &  colour )
{
    spColours[ primitive ] = colour;
}


inline void  CexmcScenePrimitives::DrawRadialLine( const G4ThreeVector &  line )
{
    radialLines.push_back( line );
    UpdateScene();
}


inline void  CexmcScenePrimitives::MarkTargetCenter( G4bool  on )
{
    markTargetCenter = on;
    UpdateScene();
}


inline void  CexmcScenePrimitives::HighlightInnerCrystals( G4bool  on )
{
    highlightInnerCrystals = on;
    UpdateScene();
}


inline void  CexmcScenePrimitives::ClearRadialLines( void )
{
    radialLines.clear();
    UpdateScene();
}


inline void CexmcScenePrimitives::UpdateScene( void )
{
    G4VVisManager *  visManager( G4VVisManager::GetConcreteInstance() );
    if ( visManager )
        visManager->NotifyHandlers();
}


#endif

