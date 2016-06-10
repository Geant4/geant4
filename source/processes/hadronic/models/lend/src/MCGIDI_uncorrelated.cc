/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include <cmath>

#include "MCGIDI_fromTOM.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
int MCGIDI_uncorrelated_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution, ptwXYPoints *norms,
        enum MCGIDI_energyType energyType, double gammaEnergy_MeV ) {

    xDataTOM_element *uncorrelatedElement;
    ptwXYPoints *angularNorms = NULL;

    if( ( energyType == MCGIDI_energyType_primaryGamma ) || ( energyType == MCGIDI_energyType_discreteGamma ) ) {
        angularNorms = norms;
        uncorrelatedElement = element; }
    else {
        if( ( uncorrelatedElement = xDataTOME_getOneElementByName( smr, element, "uncorrelated", 1 ) ) == NULL ) goto err;
    }

    if( MCGIDI_angular_parseFromTOM( smr, uncorrelatedElement, distribution, angularNorms ) ) goto err;
    if( MCGIDI_energy_parseFromTOM( smr, uncorrelatedElement, distribution, norms, energyType, gammaEnergy_MeV ) ) goto err;
    distribution->type = MCGIDI_distributionType_uncorrelated_e;

    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
int MCGIDI_uncorrelated_sampleDistribution( statusMessageReporting *smr, MCGIDI_distribution *distribution, MCGIDI_quantitiesLookupModes &modes, 
        MCGIDI_decaySamplingInfo *decaySamplingInfo ) {

    enum xDataTOM_frame frame;

    if( MCGIDI_energy_sampleEnergy( smr, distribution->energy, modes, decaySamplingInfo ) ) return( 1 );
    frame = decaySamplingInfo->frame;
    if( MCGIDI_angular_sampleMu( smr, distribution->angular, modes, decaySamplingInfo ) ) return( 1 );
    decaySamplingInfo->frame = frame;       /* Discrete and primary gammas in COM are treated as lab for now and energy sets it correctly. */
    return( 0 );
}

#if defined __cplusplus
}
#endif
