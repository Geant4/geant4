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
 *       Filename:  CexmcReconstructor.hh
 *
 *    Description:  reconstructor base class
 *
 *        Version:  1.0
 *        Created:  02.12.2009 15:44:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_RECONSTRUCTOR_HH
#define CEXMC_RECONSTRUCTOR_HH

#include <G4ThreeVector.hh>
#include <G4AffineTransform.hh>
#include "CexmcSetup.hh"
#include "CexmcCommon.hh"

class  CexmcReconstructorMessenger;
struct  CexmcEnergyDepositStore;


class  CexmcReconstructor
{
    public:
        explicit CexmcReconstructor();

        virtual ~CexmcReconstructor();

    public:
        virtual void  Reconstruct( const CexmcEnergyDepositStore *  edStore );

    public:
        void  SetCalorimeterEntryPointDefinitionAlgorithm(
                    CexmcCalorimeterEntryPointDefinitionAlgorithm  algo );

        void  SetCalorimeterEntryPointDepthDefinitionAlgorithm(
                    CexmcCalorimeterEntryPointDepthDefinitionAlgorithm  algo );

        void  SetCrystalSelectionAlgorithm(
                    CexmcCrystalSelectionAlgorithm  algo );

        void  UseInnerRefCrystal( G4bool  on = true );

        void  SetCalorimeterEntryPointDepth( G4double  depth );

        CexmcCalorimeterEntryPointDefinitionAlgorithm
                    GetCalorimeterEntryPointDefinitionAlgorithm( void ) const;

        CexmcCalorimeterEntryPointDepthDefinitionAlgorithm
                GetCalorimeterEntryPointDepthDefinitionAlgorithm( void ) const;

        CexmcCrystalSelectionAlgorithm
                    GetCrystalSelectionAlgorithm( void ) const;

        G4bool    IsInnerRefCrystalUsed( void ) const;

        G4double  GetCalorimeterEntryPointDepth( void ) const;

    public:
        const G4ThreeVector &  GetCalorimeterEPLeftPosition( void ) const;

        const G4ThreeVector &  GetCalorimeterEPRightPosition( void ) const;

        const G4ThreeVector &  GetCalorimeterEPLeftDirection( void ) const;

        const G4ThreeVector &  GetCalorimeterEPRightDirection( void ) const;

        const G4ThreeVector &  GetTargetEPPosition( void ) const;

        const G4ThreeVector &  GetTargetEPDirection( void ) const;

        const G4ThreeVector &  GetCalorimeterEPLeftWorldPosition( void ) const;

        const G4ThreeVector &  GetCalorimeterEPRightWorldPosition( void ) const;

        const G4ThreeVector &  GetCalorimeterEPLeftWorldDirection( void ) const;

        const G4ThreeVector &  GetCalorimeterEPRightWorldDirection( void )
                                                                        const;

        const G4ThreeVector &  GetTargetEPWorldPosition( void ) const;

        const G4ThreeVector &  GetTargetEPWorldDirection( void ) const;

        G4double               GetTheAngle( void ) const;

    public:
        G4bool                 HasBasicTrigger( void ) const;

        virtual G4bool         HasFullTrigger( void ) const;

    protected:
        void                   ReconstructEntryPoints(
                                    const CexmcEnergyDepositStore *  edStore );

        void                   ReconstructTargetPoint( void );

        void                   ReconstructAngle( void );

    private:
        void  CollectEDInAdjacentCrystals(
                const CexmcEnergyDepositCalorimeterCollection &  edHits,
                G4int  row, G4int  column, G4double &  ed );

        void  CalculateWeightedEPPosition(
                const CexmcEnergyDepositCalorimeterCollection &  edHits,
                G4int  row, G4int  column, G4double &  x, G4double &  y,
                G4double &  ed );

        void  TransformToAdjacentInnerCrystal( G4int &  column,
                                               G4int &  row ) const;

    protected:
        G4bool                               hasBasicTrigger;

    protected:
        CexmcCalorimeterEntryPointDefinitionAlgorithm  epDefinitionAlgorithm;

        CexmcCalorimeterEntryPointDepthDefinitionAlgorithm
                                                    epDepthDefinitionAlgorithm;

        CexmcCrystalSelectionAlgorithm       csAlgorithm;

        G4bool                               useInnerRefCrystal;

        G4double                             epDepth;

    protected:
        G4ThreeVector                        calorimeterEPLeftPosition;

        G4ThreeVector                        calorimeterEPRightPosition;

        G4ThreeVector                        calorimeterEPLeftDirection;

        G4ThreeVector                        calorimeterEPRightDirection;

        G4ThreeVector                        targetEPPosition;

        G4ThreeVector                        targetEPDirection;

        G4ThreeVector                        calorimeterEPLeftWorldPosition;

        G4ThreeVector                        calorimeterEPRightWorldPosition;

        G4ThreeVector                        calorimeterEPLeftWorldDirection;

        G4ThreeVector                        calorimeterEPRightWorldDirection;

        G4ThreeVector                        targetEPWorldPosition;

        G4ThreeVector                        targetEPWorldDirection;

        G4double                             theAngle;

    protected:
        G4double                             calorimeterEDLeftAdjacent;

        G4double                             calorimeterEDRightAdjacent;

        G4bool                               collectEDInAdjacentCrystals;

    private:
        CexmcSetup::CalorimeterGeometryData  calorimeterGeometry;

        G4AffineTransform                    calorimeterLeftTransform;
        
        G4AffineTransform                    calorimeterRightTransform;

        G4AffineTransform                    targetTransform;

        G4bool                               targetEPInitialized;

    private:
        CexmcReconstructorMessenger *        messenger;
};


inline void  CexmcReconstructor::SetCalorimeterEntryPointDefinitionAlgorithm(
                        CexmcCalorimeterEntryPointDefinitionAlgorithm  algo )
{
    epDefinitionAlgorithm = algo;
}


inline void
        CexmcReconstructor::SetCalorimeterEntryPointDepthDefinitionAlgorithm(
                    CexmcCalorimeterEntryPointDepthDefinitionAlgorithm  algo )
{
    epDepthDefinitionAlgorithm = algo;
}


inline void  CexmcReconstructor::SetCrystalSelectionAlgorithm(
                                        CexmcCrystalSelectionAlgorithm  algo )
{
    csAlgorithm = algo;
}


inline void  CexmcReconstructor::UseInnerRefCrystal( G4bool  on )
{
    useInnerRefCrystal = on;
}


inline void  CexmcReconstructor::SetCalorimeterEntryPointDepth(
                                                            G4double  depth )
{
    epDepth = depth;
}


inline CexmcCalorimeterEntryPointDefinitionAlgorithm
        CexmcReconstructor::GetCalorimeterEntryPointDefinitionAlgorithm( void )
                                                                        const
{
    return epDefinitionAlgorithm;
}


inline CexmcCalorimeterEntryPointDepthDefinitionAlgorithm
        CexmcReconstructor::GetCalorimeterEntryPointDepthDefinitionAlgorithm(
                                                                    void ) const
{
    return epDepthDefinitionAlgorithm;
}


inline CexmcCrystalSelectionAlgorithm
                CexmcReconstructor::GetCrystalSelectionAlgorithm( void ) const
{
    return csAlgorithm;
}


inline G4bool  CexmcReconstructor::IsInnerRefCrystalUsed( void ) const
{
    return useInnerRefCrystal;
}


inline G4double  CexmcReconstructor::GetCalorimeterEntryPointDepth( void ) const
{
    return epDepth;
}


inline const G4ThreeVector &
                CexmcReconstructor::GetCalorimeterEPLeftPosition( void ) const
{
    return calorimeterEPLeftPosition;
}


inline const G4ThreeVector &
                CexmcReconstructor::GetCalorimeterEPRightPosition( void ) const
{
    return calorimeterEPRightPosition;
}


inline const G4ThreeVector &
                CexmcReconstructor::GetCalorimeterEPLeftDirection( void ) const
{
    return calorimeterEPLeftDirection;
}


inline const G4ThreeVector &
                CexmcReconstructor::GetCalorimeterEPRightDirection( void ) const
{
    return calorimeterEPRightDirection;
}


inline const G4ThreeVector &
                CexmcReconstructor::GetTargetEPPosition( void ) const
{
    return targetEPPosition;
}


inline const G4ThreeVector &
                CexmcReconstructor::GetTargetEPDirection( void ) const
{
    return targetEPDirection;
}


inline const G4ThreeVector &
        CexmcReconstructor::GetCalorimeterEPLeftWorldPosition( void ) const
{
    return calorimeterEPLeftWorldPosition;
}


inline const G4ThreeVector &
        CexmcReconstructor::GetCalorimeterEPRightWorldPosition( void ) const
{
    return calorimeterEPRightWorldPosition;
}


inline const G4ThreeVector &
        CexmcReconstructor::GetCalorimeterEPLeftWorldDirection( void ) const
{
    return calorimeterEPLeftWorldDirection;
}


inline const G4ThreeVector &
        CexmcReconstructor::GetCalorimeterEPRightWorldDirection( void ) const
{
    return calorimeterEPRightWorldDirection;
}


inline const G4ThreeVector &
        CexmcReconstructor::GetTargetEPWorldPosition( void ) const
{
    return targetEPWorldPosition;
}


inline const G4ThreeVector &
        CexmcReconstructor::GetTargetEPWorldDirection( void ) const
{
    return targetEPWorldDirection;
}


inline G4double  CexmcReconstructor::GetTheAngle( void ) const
{
    return theAngle;
}


inline G4bool  CexmcReconstructor::HasBasicTrigger( void ) const
{
    return hasBasicTrigger;
}


inline void  CexmcReconstructor::TransformToAdjacentInnerCrystal(
                                        G4int &  column, G4int &  row ) const
{
    if ( column == 0 )
        ++column;
    if ( column == calorimeterGeometry.nCrystalsInRow - 1 )
        --column;
    if ( row == 0 )
        ++row;
    if ( row == calorimeterGeometry.nCrystalsInColumn - 1 )
        --row;
}


#endif

