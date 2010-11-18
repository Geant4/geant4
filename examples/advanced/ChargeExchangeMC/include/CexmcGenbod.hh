/*
 * =============================================================================
 *
 *       Filename:  CexmcGenbod.hh
 *
 *    Description:  original fortran routine GENBOD wrapper
 *
 *        Version:  1.0
 *        Created:  08.09.2010 14:32:14
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_GENBOD_HH
#define CEXMC_GENBOD_HH

#ifdef CEXMC_USE_GENBOD

#include "CexmcPhaseSpaceGenerator.hh"


class  CexmcGenbod : public CexmcPhaseSpaceGenerator
{
    public:
        CexmcGenbod();

    public:
        G4bool    CheckKinematics( void );

        G4double  Generate( void );

    private:
        void      ParticleChangeHook( void );

        void      FermiEnergyDepStatusChangeHook( void );
};

#endif

#endif

