/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include <cmath>

#include "MCGIDI.h"
#include "MCGIDI_misc.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
MCGIDI_distribution *MCGIDI_distribution_new( statusMessageReporting *smr ) {

    MCGIDI_distribution *distribution;

    if( ( distribution = (MCGIDI_distribution *) smr_malloc2( smr, sizeof( MCGIDI_distribution ), 0, "distribution" ) ) == NULL ) return( NULL );
    if( MCGIDI_distribution_initialize( smr, distribution ) ) distribution = MCGIDI_distribution_free( smr, distribution );
    return( distribution );
}
/*
************************************************************
*/
int MCGIDI_distribution_initialize( statusMessageReporting * /*smr*/, MCGIDI_distribution *distribution ) {

    memset( distribution, 0, sizeof( MCGIDI_distribution ) );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_distribution *MCGIDI_distribution_free( statusMessageReporting *smr, MCGIDI_distribution *distribution ) {

    MCGIDI_distribution_release( smr, distribution );
    smr_freeMemory( (void **) &distribution );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_distribution_release( statusMessageReporting *smr, MCGIDI_distribution *distribution ) {

    if( distribution->angular ) distribution->angular = MCGIDI_angular_free( smr, distribution->angular );
    if( distribution->energy ) distribution->energy = MCGIDI_energy_free( smr, distribution->energy );
    if( distribution->KalbachMann ) distribution->KalbachMann = MCGIDI_KalbachMann_free( smr, distribution->KalbachMann );    
    if( distribution->energyAngular ) distribution->energyAngular = MCGIDI_energyAngular_free( smr, distribution->energyAngular );
    if( distribution->angularEnergy ) distribution->angularEnergy = MCGIDI_angularEnergy_free( smr, distribution->angularEnergy );

    MCGIDI_distribution_initialize( smr, distribution );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_distribution_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_product *product, MCGIDI_POPs * /*pops*/, ptwXYPoints *norms ) {

    char const *nativeData, *gammaEnergy;
    double gammaEnergy_MeV;
    MCGIDI_distribution *distribution = &(product->distribution);
    xDataTOM_element *distributionElement;
    enum MCGIDI_energyType energyType = MCGIDI_energyType_unknown;

    MCGIDI_distribution_initialize( smr, distribution );
 
    distribution->product = product;
    if( ( distributionElement = xDataTOME_getOneElementByName( smr, element, "distributions", 1 ) ) == NULL ) goto err;
    if( ( nativeData = xDataTOM_getAttributesValueInElement( distributionElement, "nativeData" ) ) == NULL ) goto err;

    if( strcmp( product->pop->name, "gamma" ) == 0 ) {
        if( ( gammaEnergy = xDataTOM_getAttributesValueInElement( element, "discrete" ) ) != NULL ) {
            if( MCGIDI_misc_PQUStringToDoubleInUnitOf( smr, gammaEnergy, "MeV", &gammaEnergy_MeV ) ) goto err;
            energyType = MCGIDI_energyType_discreteGamma; }
        else if( ( gammaEnergy = xDataTOM_getAttributesValueInElement( element, "primary" ) ) != NULL ) {
            if( MCGIDI_misc_PQUStringToDoubleInUnitOf( smr, gammaEnergy, "MeV", &gammaEnergy_MeV ) ) goto err;
            energyType = MCGIDI_energyType_primaryGamma;
        }
        if( gammaEnergy != NULL ) {
            if( strcmp( nativeData, "angular" ) ) {
                smr_setReportError2( smr, smr_unknownID, 1, "%s gamma can only have a distribution with 'nativeData' = 'angular' and not '%s'",
                    gammaEnergy, nativeData );
                goto err;
            }
            nativeData = "uncorrelated";
        }
    }

    if( strcmp( nativeData, "angular" ) == 0 ) {
        if( MCGIDI_angular_parseFromTOM( smr, distributionElement, distribution, norms ) ) goto err; }
    else if( strcmp( nativeData, "uncorrelated" ) == 0 ) {
        if( MCGIDI_uncorrelated_parseFromTOM( smr, distributionElement, distribution, norms, energyType, gammaEnergy_MeV ) ) goto err; }
    else if( strcmp( nativeData, "energyAngular" ) == 0 ) {
        if( MCGIDI_energyAngular_parseFromTOM( smr, distributionElement, distribution ) ) goto err; }
    else if( strcmp( nativeData, "angularEnergy" ) == 0 ) {
        if( MCGIDI_angularEnergy_parseFromTOM( smr, distributionElement, distribution ) ) goto err; }
    else if( strcmp( nativeData, "Legendre" ) == 0 ) {
        if( MCGIDI_energyAngular_parseFromTOM( smr, distributionElement, distribution ) ) goto err; }
    else if( strcmp( nativeData, "LLNLAngular_angularEnergy" ) == 0 ) {
        if( MCGIDI_LLNLAngular_angularEnergy_parseFromTOM( smr, distributionElement, distribution ) ) goto err; }
    else if( strcmp( nativeData, "none" ) == 0 ) {
        distribution->type = MCGIDI_distributionType_none_e; }
    else if( strcmp( nativeData, "unknown" ) == 0 ) {
        distribution->type = MCGIDI_distributionType_unknown_e; }
    else {
        smr_setReportError2( smr, smr_unknownID, 1, "Unsupported distribution = '%s'\n", nativeData );
        goto err;
    }

    return( 0 );

err:
    MCGIDI_distribution_release( smr, distribution );
    return( 1 );
}

#if defined __cplusplus
}
#endif

