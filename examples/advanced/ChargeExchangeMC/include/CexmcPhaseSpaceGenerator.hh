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
 *       Filename:  CexmcPhaseSpaceGenerator.hh
 *
 *    Description:  phase space generator interface
 *
 *        Version:  1.0
 *        Created:  08.09.2010 13:42:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PHASE_SPACE_GENERATOR_HH
#define CEXMC_PHASE_SPACE_GENERATOR_HH

#include <vector>
#include <G4Types.hh>
#include <G4LorentzVector.hh>


struct  CexmcPhaseSpaceOutVectorElement
{
    CexmcPhaseSpaceOutVectorElement( G4LorentzVector *  lVec_,
                                     G4double  mass_ ) :
        lVec( lVec_ ), mass( mass_ )
    {}

    G4LorentzVector *  lVec;

    G4double           mass;
};


typedef std::vector< const G4LorentzVector * >  CexmcPhaseSpaceInVector;

typedef std::vector< CexmcPhaseSpaceOutVectorElement >
                                                CexmcPhaseSpaceOutVector;


class  CexmcPhaseSpaceGenerator
{
    public:
        CexmcPhaseSpaceGenerator();

        virtual ~CexmcPhaseSpaceGenerator();

    public:
        virtual G4bool    CheckKinematics( void );

        virtual G4double  Generate( void ) = 0;

    public:
        void  SetParticles( const CexmcPhaseSpaceInVector &  inVec_,
                            const CexmcPhaseSpaceOutVector &  outVec_ );

        void  SetFermiEnergyDependence( G4bool  on = true );

    protected:
        virtual void    ParticleChangeHook( void );

        virtual void    FermiEnergyDepStatusChangeHook( void );

    protected:
        CexmcPhaseSpaceInVector   inVec;

        CexmcPhaseSpaceOutVector  outVec;

        G4bool                    fermiEnergyDepIsOn;

        G4double                  totalEnergy;

        G4double                  totalMass;
};


#endif

