/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include <cmath>

#include "MCGIDI.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
int MCGIDI_kinetics_2BodyReaction( statusMessageReporting *smr, MCGIDI_angular *angular, double K, double mu, double phi, 
        MCGIDI_sampledProductsData *outgoingData ) {

    double m1 = angular->projectileMass_MeV, m2 = angular->targetMass_MeV, m3 = angular->productMass_MeV, m4 = angular->residualMass_MeV, mi, mf, Kp, x, beta;

    mi = m1 + m2;
    mf = m3 + m4;
    beta = std::sqrt( K * ( K + 2. * m1 ) ) / ( K + mi );
    x = K * m2 / ( mi * mi );
    if( x < 2e-5 ) {                                        /* Kp is the total kinetic energy for m3 and m4 in the COM frame. */
        Kp = mi - mf + K * m2 / mi * ( 1 - 0.5 * x * ( 1 - x ) ); }
    else {
        Kp = std::sqrt( mi * mi + 2 * K * m2 ) - mf;
    }
    if( Kp < 0 ) Kp = 0.;           /* ???? There needs to be a better test here. */
    return( MCGIDI_kinetics_COMKineticEnergy2LabEnergyAndMomentum( smr, beta, Kp, mu, phi, m3, m4, outgoingData ) );
}
/*
************************************************************
*/
int MCGIDI_kinetics_COMKineticEnergy2LabEnergyAndMomentum( statusMessageReporting * /*smr*/, double beta, double e_kinetic_com, double mu, double phi, 
        double m3cc, double m4cc, MCGIDI_sampledProductsData *outgoingData ) {
/*
*   beta            the velocity/speedOflight of the com frame relative to the lab frame.
*   e_kinetic_com   Total kinetic energy (K1 + K2) in the COM frame.
*   mu              cos( theta ) in the COM frame.
*/
    double x, v_p, p, pp3, pp4, px3, py3, pz3, pz4, pz, p_perp2, E3, E4, gamma, m3cc2 = m3cc * m3cc, m4cc2 = m4cc * m4cc;

    p = std::sqrt( e_kinetic_com * ( e_kinetic_com + 2. * m3cc ) * ( e_kinetic_com + 2. * m4cc )  * ( e_kinetic_com + 2. * ( m3cc + m4cc ) ) ) /
            ( 2. * ( e_kinetic_com + m3cc + m4cc ) );
    py3 = p * std::sqrt( 1 - mu * mu );
    px3 = py3 * std::cos( phi );
    py3 *= std::sin( phi );
    pz = p * mu;
    if( 1 ) {                           /* ????????? Assuming the answer is wanted in the lab frame for now. */
        E3 = std::sqrt( p * p + m3cc2 );
        E4 = std::sqrt( p * p + m4cc2 );
        gamma = std::sqrt( 1. / ( 1. - beta * beta ) );
        pz3 = gamma * (  pz + beta * E3 );
        pz4 = gamma * ( -pz + beta * E4 ); }
    else {                              /* COM frame. */
        pz3 = pz;
        pz4 = -pz;
    }
    outgoingData[1].isVelocity = outgoingData[0].isVelocity;

    p_perp2 = px3 * px3 + py3 * py3;

    outgoingData[0].px_vx = px3;
    outgoingData[0].py_vy = py3;
    outgoingData[0].pz_vz = pz3;
    pp3 = p_perp2 + pz3 * pz3;
    x = ( m3cc > 0 ) ? pp3 / ( 2 * m3cc2 ) : 1.;
    if( x < 1e-5 ) {
        outgoingData[0].kineticEnergy = m3cc * x  * ( 1 - 0.5 * x * ( 1 - x ) ); }
    else {
        outgoingData[0].kineticEnergy = std::sqrt( m3cc2 + pp3 ) - m3cc;
    }
    outgoingData[1].px_vx = -px3;
    outgoingData[1].py_vy = -py3;
    outgoingData[1].pz_vz = pz4;
    pp4 = p_perp2 + pz4 * pz4;
    x = ( m4cc > 0 ) ? pp4 / ( 2 * m4cc2 ) : 1.;
    if( x < 1e-5 ) {
        outgoingData[1].kineticEnergy = m4cc * x  * ( 1 - 0.5 * x * ( 1 - x ) ); }
    else {
        outgoingData[1].kineticEnergy = std::sqrt( m4cc2 + pp4 ) - m4cc;
    }

    if( outgoingData[0].isVelocity ) {
        v_p = MCGIDI_speedOfLight_cm_sec / std::sqrt( pp3 + m3cc2 );
        outgoingData[0].px_vx *= v_p;
        outgoingData[0].py_vy *= v_p;
        outgoingData[0].pz_vz *= v_p;

        v_p = MCGIDI_speedOfLight_cm_sec / std::sqrt( pp4 + m4cc2 );
        outgoingData[1].px_vx *= v_p;
        outgoingData[1].py_vy *= v_p;
        outgoingData[1].pz_vz *= v_p;
    }

    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_kinetics_COM2Lab( statusMessageReporting *smr, MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo, double masses[3] ) {
/*      
*       massProjectile = masses[0], massTarget = masses[1], massProduct = masses[2];
*/
    double a = masses[0] + masses[1], b, e_in = modes.getProjectileEnergy( ) * masses[0] * masses[2] / ( a * a ), Ep;

    if( decaySamplingInfo->frame != xDataTOM_frame_centerOfMass ) {
        smr_setReportError2( smr, smr_unknownID, 1, "bad frame = %d for COM to lab conversion of mu/energy", decaySamplingInfo->frame );
        return( 1 );
    }
    a = std::sqrt( e_in );
    b = std::sqrt( decaySamplingInfo->Ep );
    Ep = decaySamplingInfo->Ep + e_in + 2. * decaySamplingInfo->mu * a * b;
    if( Ep != 0 ) {
        decaySamplingInfo->mu = ( a + decaySamplingInfo->mu * b ) / std::sqrt( Ep );
    }
    decaySamplingInfo->Ep = Ep;
    decaySamplingInfo->frame = xDataTOM_frame_lab;
    return( 0 );
}

#if defined __cplusplus
}
#endif
