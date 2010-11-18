/*
 * =============================================================================
 *
 *       Filename:  CexmcGenbod.cc
 *
 *    Description:  original fortran routine GENBOD wrapper
 *
 *        Version:  1.0
 *        Created:  08.09.2010 14:36:23
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifdef CEXMC_USE_GENBOD

#include "CexmcGenbod.hh"
#include "CexmcException.hh"


extern "C"
{
    extern int  genbod_( void );
}


extern struct  genbod_in_data
{
    int    np;
    float  tecm;
    float  amass[ 18 ];
    int    kgenev;
}  genin_;


extern struct  genbod_out_data
{
    float  pcm[ 18 ][ 5 ];
    float  wt;
}  genout_;


CexmcGenbod::CexmcGenbod()
{
    genin_.np = 0;
    genin_.kgenev = 1;
}


G4bool  CexmcGenbod::CheckKinematics( void )
{
    totalEnergy = 0.;
    for ( CexmcPhaseSpaceInVector::const_iterator  k( inVec.begin() );
                                                        k != inVec.end(); ++k )
    {
        totalEnergy += ( *k )->e();
    }

    /* epsilon is needed to compensate float to fortran real cast accuracy,
     * the value 3E-6 was found experimentally, maybe has to be made bigger,
     * but not smaller */
    const float  epsilon( 3E-6 );

    return totalEnergy - totalMass > 0.0f + epsilon;
}


G4double  CexmcGenbod::Generate( void )
{
    genin_.tecm = totalEnergy / GeV;

    genbod_();

    float ( *pcm )[ 5 ]( &genout_.pcm[ 0 ] );

    for ( CexmcPhaseSpaceOutVector::iterator  k( outVec.begin() );
                                                        k != outVec.end(); ++k )
    {
        k->lVec->setPx( ( *pcm )[ 0 ] * GeV );
        k->lVec->setPy( ( *pcm )[ 1 ] * GeV );
        k->lVec->setPz( ( *pcm )[ 2 ] * GeV );
        k->lVec->setE( ( *pcm++ )[ 3 ] * GeV );
    }

    return genout_.wt;
}


void  CexmcGenbod::ParticleChangeHook( void )
{
    size_t  nmbOfOutputParticles( outVec.size() );

    if ( nmbOfOutputParticles < 2 || nmbOfOutputParticles > 18 )
        throw CexmcException( CexmcKinematicsException );

    genin_.np = nmbOfOutputParticles;

    float *  amass( genin_.amass );

    for ( CexmcPhaseSpaceOutVector::const_iterator  k( outVec.begin() );
                                                        k != outVec.end(); ++k )
    {
        *amass++ = k->mass / GeV;
    }
}


void  CexmcGenbod::FermiEnergyDepStatusChangeHook( void )
{
    genin_.kgenev = fermiEnergyDepIsOn ? 2 : 1;
}

#endif

