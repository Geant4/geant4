/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include "MCGIDI.h"

using namespace GIDI;
/*  ---- MCGIDI_samplingMethods ----  */
/*
=========================================================
*/
MCGIDI_samplingMethods::MCGIDI_samplingMethods( ) {

}
/*
=========================================================
*/
MCGIDI_samplingMethods::~MCGIDI_samplingMethods( ) {

}

/*  ---- MCGIDI_samplingSettings ----  */
/*
=========================================================
*/
MCGIDI_samplingSettings::MCGIDI_samplingSettings( enum xDataTOM_frame frame, bool wantVelocities, double (*rng)( void * ), void *rngState ) {

    mWantFrame = frame;
    mWantVelocities = wantVelocities;
    mRng = rng;
    mRngState = rngState;

    mGotFrame = xDataTOM_frame_invalid;
    mPoP = NULL;
    mMu = 0.;
    mEp = 0.; 
}
/*
=========================================================
*/
MCGIDI_samplingSettings::~MCGIDI_samplingSettings( void ) {

}
/*
=========================================================
*/
int MCGIDI_samplingSettings::setProductMultiplicityBias( statusMessageReporting *smr, int PoPID, double factor ) {

    if( factor < 0 ) {
        smr_setReportError2( smr, smr_unknownID, 1, "factor = %e cannot be negative", factor );
        return( 1 );
    }

    for( int i1 = 0; i1 < (int) mSamplingMultiplicityBiases.size( ); ++i1 ) {
        if( PoPID == mSamplingMultiplicityBiases[i1].PoPID ) {
            mSamplingMultiplicityBiases[i1].multiplicityFactor = factor;
            return( 0 );
        }
    }
    MCGIDI_samplingMultiplicityBias samplingMultiplicityBias = { PoPID, factor };
    mSamplingMultiplicityBiases.push_back( samplingMultiplicityBias );
    return( 0 );
}
