/*
 * =============================================================================
 *
 *       Filename:  CexmcReimplementedGenbod.hh
 *
 *    Description:  reimplemented GENBOD
 *                  (mostly adopted from ROOT TGenPhaseSpace)
 *
 *        Version:  1.0
 *        Created:  08.09.2010 18:46:18
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_REIMPLEMENTED_GENBOD_HH
#define CEXMC_REIMPLEMENTED_GENBOD_HH

#include "CexmcPhaseSpaceGenerator.hh"


class  CexmcReimplementedGenbod : public CexmcPhaseSpaceGenerator
{
    public:
        CexmcReimplementedGenbod();

    public:
        G4double  Generate( void );

    private:
        void      ParticleChangeHook( void );

        void      FermiEnergyDepStatusChangeHook( void );

    private:
        void      SetMaxWeight( void );

        G4double  PDK( G4double  a, G4double  b, G4double  c );

    private:
        G4double  maxWeight;

        G4int     nmbOfOutputParticles;

    private:
        static const G4int  maxParticles = 18;
};


#endif

