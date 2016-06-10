/*
# <<BEGIN-copyright>>
# Copyright (c) 2010, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory 
# Written by Bret R. Beck, beck6@llnl.gov. 
# CODE-461393
# All rights reserved. 
#  
# This file is part of GIDI. For details, see nuclear.llnl.gov. 
# Please also read the "Additional BSD Notice" at nuclear.llnl.gov. 
# 
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met: 
#
#      1) Redistributions of source code must retain the above copyright notice, 
#         this list of conditions and the disclaimer below.
#      2) Redistributions in binary form must reproduce the above copyright notice, 
#         this list of conditions and the disclaimer (as noted below) in the 
#          documentation and/or other materials provided with the distribution.
#      3) Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#         used to endorse or promote products derived from this software without 
#         specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT 
# SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS 
# OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# <<END-copyright>>
*/
#include <iostream>
#include <float.h>
#include <string.h>
#include <cmath>
#include <tpia_target.h>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif 

/*
************************************************************
*/
int tpia_kinetics_2BodyReaction( statusMessageReporting *smr, tpia_decayChannel *decayChannel, double K, double mu, double phi, 
        tpia_productOutgoingData *outgoingData ) {

    tpia_product *pp3 = tpia_decayChannel_getFirstProduct( decayChannel ), *pp4;
    double m1 = decayChannel->m1_fullMass_MeV, m2 = decayChannel->m2_fullMass_MeV, m3, m4, mi, mf, Kp, x, beta;

    pp4 = tpia_decayChannel_getNextProduct( pp3 );
    m3 = pp3->productID->fullMass_MeV;
    m4 = pp4->productID->fullMass_MeV;
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
    outgoingData[0].decayChannel = &(pp3->decayChannel);
    outgoingData[1].genre = outgoingData[0].genre;
    outgoingData[1].productID = pp4->productID;
    outgoingData[1].decayChannel = &(pp4->decayChannel);
    return( tpia_kinetics_COMKineticEnergy2LabEnergyAndMomentum( smr, beta, Kp, mu, phi, m3, m4, outgoingData ) );
}
/*
************************************************************
*/
//int tpia_kinetics_COMKineticEnergy2LabEnergyAndMomentum( statusMessageReporting *smr, double beta, double e_kinetic_com, double mu, double phi, 
int tpia_kinetics_COMKineticEnergy2LabEnergyAndMomentum( statusMessageReporting *, double beta, double e_kinetic_com, double mu, double phi, 
        double m3cc, double m4cc, tpia_productOutgoingData *outgoingData ) {
/*
*   beta            the velocity/speedOflight of the com frame relative to the lab frame.
*   e_kinetic_com   Total kinetic energy (K1 + K2) in the COM frame.
*   mu              std::cos( theta ) in the COM frame.
*/
    double x, v_p, p, pp3, pp4, px3, py3, pz3, pz4, pz, p_perp2, E3, E4, gamma, m3cc2 = m3cc * m3cc, m4cc2 = m4cc * m4cc;

    p = std::sqrt( e_kinetic_com * ( e_kinetic_com + 2. * m3cc ) * ( e_kinetic_com + 2. * m4cc )  * ( e_kinetic_com + 2. * ( m3cc + m4cc ) ) ) /
            ( 2. * ( e_kinetic_com + m3cc + m4cc ) );
    py3 = p * std::sqrt( 1 - mu * mu );
    px3 = py3 * std::cos( phi );
    py3 *= std::sin( phi );
    pz = p * mu;
    if( tpia_frame_getColumn( NULL, &(outgoingData[0].frame), 0 ) == tpia_referenceFrame_lab ) {
        E3 = std::sqrt( p * p + m3cc2 );
        E4 = std::sqrt( p * p + m4cc2 );
        gamma = std::sqrt( 1. / ( 1. - beta * beta ) );
        pz3 = gamma * (  pz + beta * E3 );
        pz4 = gamma * ( -pz + beta * E4 ); }
    else {
        pz3 = pz;
        pz4 = -pz;
    }
    outgoingData[1].isVelocity = outgoingData[0].isVelocity;
    outgoingData[1].frame = outgoingData[0].frame;

    p_perp2 = px3 * px3 + py3 * py3;

    outgoingData[0].px_vx = px3;
    outgoingData[0].py_vy = py3;
    outgoingData[0].pz_vz = pz3;
    pp3 = p_perp2 + pz3 * pz3;
//TK140602 Modified for protecting divided by 0 BEGIN
    if ( m3cc2 != 0 ) 
    x = pp3 / ( 2 * m3cc2 );
    else
    x = FLT_MIN;
//TK140602 Modified for protecting divided by 0 END
    if( x < 1e-5 ) {
        outgoingData[0].kineticEnergy = m3cc * x  * ( 1 - 0.5 * x * ( 1 - x ) ); }
    else {
        outgoingData[0].kineticEnergy = std::sqrt( m3cc2 + pp3 ) - m3cc;
    }
    outgoingData[1].px_vx = -px3;
    outgoingData[1].py_vy = -py3;
    outgoingData[1].pz_vz = pz4;
    pp4 = p_perp2 + pz4 * pz4;
    x = pp4 / ( 2 * m4cc2 );
    if( x < 1e-5 ) {
        outgoingData[1].kineticEnergy = m4cc * x  * ( 1 - 0.5 * x * ( 1 - x ) ); }
    else {
        outgoingData[1].kineticEnergy = std::sqrt( m4cc2 + pp4 ) - m4cc;
    }

    if( outgoingData[0].isVelocity ) {
        v_p = tpia_speedOfLight_cm_sec / std::sqrt( pp3 + m3cc2 );
        outgoingData[0].px_vx *= v_p;
        outgoingData[0].py_vy *= v_p;
        outgoingData[0].pz_vz *= v_p;

        v_p = tpia_speedOfLight_cm_sec / std::sqrt( pp4 + m4cc2 );
        outgoingData[1].px_vx *= v_p;
        outgoingData[1].py_vy *= v_p;
        outgoingData[1].pz_vz *= v_p;
    }

    return( 0 );
}

#if defined __cplusplus
}
#endif
