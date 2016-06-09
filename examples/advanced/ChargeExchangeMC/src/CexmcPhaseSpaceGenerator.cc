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
 *       Filename:  CexmcPhaseSpaceGenerator.cc
 *
 *    Description:  phase space generator interface
 *
 *        Version:  1.0
 *        Created:  08.09.2010 14:05:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include "CexmcPhaseSpaceGenerator.hh"


CexmcPhaseSpaceGenerator::CexmcPhaseSpaceGenerator() :
    fermiEnergyDepIsOn( false ), totalEnergy( 0. ), totalMass( 0. )
{
}


CexmcPhaseSpaceGenerator::~CexmcPhaseSpaceGenerator()
{
}


void  CexmcPhaseSpaceGenerator::SetParticles(
                                    const CexmcPhaseSpaceInVector &  inVec_,
                                    const CexmcPhaseSpaceOutVector &  outVec_ )
{
    inVec = inVec_;
    outVec = outVec_;
    for ( CexmcPhaseSpaceOutVector::const_iterator  k( outVec.begin() );
                                                        k != outVec.end(); ++k )
    {
        totalMass += k->mass;
    }

    ParticleChangeHook();
}


void  CexmcPhaseSpaceGenerator::SetFermiEnergyDependence( G4bool  on )
{
    fermiEnergyDepIsOn = on;

    FermiEnergyDepStatusChangeHook();
}


G4bool  CexmcPhaseSpaceGenerator::CheckKinematics( void )
{
    totalEnergy = 0.;
    for ( CexmcPhaseSpaceInVector::const_iterator  k( inVec.begin() );
                                                        k != inVec.end(); ++k )
    {
        totalEnergy += ( *k )->e();
    }

    return totalEnergy > totalMass;
}


void  CexmcPhaseSpaceGenerator::ParticleChangeHook( void )
{
}


void  CexmcPhaseSpaceGenerator::FermiEnergyDepStatusChangeHook( void )
{
}

