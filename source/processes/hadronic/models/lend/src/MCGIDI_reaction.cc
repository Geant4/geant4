/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include <cmath>

#include <PoPs.h>
#include "MCGIDI.h"
#include "MCGIDI_misc.h"
#include "MCGIDI_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

#define nParticleChanges 6

static int MCGIDI_reaction_initialize2( statusMessageReporting *smr, MCGIDI_reaction *reaction );
static int MCGIDI_reaction_particleChanges( MCGIDI_POP *projectile, MCGIDI_POP *target, MCGIDI_productsInfo *productsInfo, int n1, int *particlesChanges );
static int MCGIDI_reaction_ParseReactionTypeAndDetermineProducts( statusMessageReporting *smr, MCGIDI_POPs *pops, MCGIDI_reaction *reaction );
static int MCGIDI_reaction_ParseDetermineReactionProducts( statusMessageReporting *smr, MCGIDI_POPs *pops, MCGIDI_outputChannel *outputChannel,
    MCGIDI_productsInfo *productsInfo, MCGIDI_reaction *reaction, double *finalQ, int level );
static int MCGIDI_reaction_addReturnProduct( statusMessageReporting *smr, MCGIDI_productsInfo *productsInfo, int ID, MCGIDI_product *product, 
        MCGIDI_reaction *reaction, int transportable );
static int MCGIDI_reaction_setENDL_CSNumbers( statusMessageReporting *smr, MCGIDI_reaction *reaction );
/*
************************************************************
*/
MCGIDI_reaction *MCGIDI_reaction_new( statusMessageReporting *smr ) {

    MCGIDI_reaction *reaction;

    if( ( reaction = (MCGIDI_reaction *) smr_malloc2( smr, sizeof( MCGIDI_reaction ), 0, "reaction" ) ) == NULL ) return( NULL );
    if( MCGIDI_reaction_initialize( smr, reaction ) ) reaction = MCGIDI_reaction_free( smr, reaction );
    return( reaction );
}
/*
************************************************************
*/
int MCGIDI_reaction_initialize( statusMessageReporting *smr, MCGIDI_reaction *reaction ) {

    if( MCGIDI_reaction_initialize2( smr, reaction ) != 0 ) return( 1 );
    reaction->transportabilities = new transportabilitiesMap( );
    return( 0 );
}
/*
************************************************************
*/
static int MCGIDI_reaction_initialize2( statusMessageReporting *smr, MCGIDI_reaction *reaction ) {

    memset( reaction, 0, sizeof( MCGIDI_reaction ) );
    xDataTOMAL_initial( smr, &(reaction->attributes) );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_reaction *MCGIDI_reaction_free( statusMessageReporting *smr, MCGIDI_reaction *reaction ) {

    MCGIDI_reaction_release( smr, reaction );
    smr_freeMemory( (void **) &reaction );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_reaction_release( statusMessageReporting *smr, MCGIDI_reaction *reaction ) {

    ptwXY_free( reaction->crossSection );
    ptwX_free( reaction->crossSectionGrouped );
    MCGIDI_outputChannel_release( smr, &(reaction->outputChannel) );
    xDataTOMAL_release( &(reaction->attributes) );
    smr_freeMemory( (void **) &(reaction->outputChannelStr) );
    if( reaction->productsInfo.productInfo != NULL ) smr_freeMemory( (void **) &(reaction->productsInfo.productInfo) );
    delete reaction->transportabilities;
    MCGIDI_reaction_initialize2( smr, reaction );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_reaction_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_target_heated *target, 
        MCGIDI_POPs *pops, MCGIDI_reaction *reaction ) {

    xDataTOM_element *child, *linear, *outputChannel;
    enum xDataTOM_interpolationFlag independent, dependent;
    enum xDataTOM_interpolationQualifier qualifier;
    char const *outputChannelStr, *crossSectionUnits[2] = { "MeV", "b" };

    MCGIDI_reaction_initialize( smr, reaction );

    reaction->target = target;
    reaction->reactionType = MCGIDI_reactionType_unknown_e;
    if( xDataTOME_copyAttributionList( smr, &(reaction->attributes), element ) ) goto err;
    if( xDataTOME_convertAttributeToInteger( smr, element, "ENDF_MT", &(reaction->ENDF_MT) ) ) goto err;
    if( ( outputChannelStr = xDataTOM_getAttributesValueInElement( element, "outputChannel" ) ) == NULL ) goto err;
    if( ( reaction->outputChannelStr = smr_allocateCopyString2( smr, outputChannelStr, "reaction->outputChannelStr" ) ) == NULL ) goto err;

    if( ( child = xDataTOME_getOneElementByName( smr, element, "crossSection", 1 ) ) == NULL ) goto err;
    if( ( linear = xDataTOME_getOneElementByName( smr, child, "linear", 0 ) ) == NULL ) {
        if( ( linear = xDataTOME_getOneElementByName( smr, child, "pointwise", 1 ) ) == NULL ) goto err;
    }
    if( xDataTOME_getInterpolation( smr, linear, 0, &independent, &dependent, &qualifier ) ) goto err;
    if( ( independent != xDataTOM_interpolationFlag_linear ) || ( dependent != xDataTOM_interpolationFlag_linear ) ) {
        smr_setReportError2( smr, smr_unknownID, 1, "cross section interpolation (%d,%d) is not linear-linear", independent, dependent );
        goto err;
    }
    if( ( reaction->crossSection = MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( smr, linear, crossSectionUnits ) ) == NULL ) goto err;
    reaction->domainValuesPresent = 1;
    reaction->EMin = ptwXY_getXMin( reaction->crossSection );
    reaction->EMax = ptwXY_getXMax( reaction->crossSection );

    if( ( outputChannel = xDataTOME_getOneElementByName( smr, element, "outputChannel", 1 ) ) == NULL ) goto err;
    if( MCGIDI_outputChannel_parseFromTOM( smr, outputChannel, pops, &(reaction->outputChannel), reaction, NULL ) ) goto err;

    if( MCGIDI_reaction_ParseReactionTypeAndDetermineProducts( smr, pops, reaction ) != 0 ) goto err;

    return( 0 );

err:
    MCGIDI_reaction_release( smr, reaction );
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_reaction_ParseReactionTypeAndDetermineProducts( statusMessageReporting *smr, MCGIDI_POPs *pops, MCGIDI_reaction *reaction ) {

    MCGIDI_outputChannel *outputChannel = &(reaction->outputChannel);
    int MT;
    int particlesChanges[nParticleChanges], numberOfChanges;
    double finalQ = 0.;

    if( MCGIDI_reaction_ParseDetermineReactionProducts( smr, pops, outputChannel, &(reaction->productsInfo), reaction, &finalQ, 0 ) != 0 ) return( 1 );
    reaction->finalQ = finalQ;
    MT = MCGIDI_reaction_getENDF_MTNumber( reaction );
    switch( MT ) {
    case 2 :
        reaction->reactionType = MCGIDI_reactionType_elastic_e;
        break;
    case 18 : case 19 : case 20 : case 21 : case 38 :
        reaction->reactionType = MCGIDI_reactionType_fission_e;
        break;
    case 102 :
        reaction->reactionType = MCGIDI_reactionType_capture_e;
        break;
    case 5 :
        reaction->reactionType = MCGIDI_reactionType_sumOfRemainingOutputChannels_e;
        break;
    default :
        numberOfChanges = MCGIDI_reaction_particleChanges( reaction->target->projectilePOP, reaction->target->targetPOP, &(reaction->productsInfo), 
            nParticleChanges, particlesChanges );

        reaction->reactionType = MCGIDI_reactionType_unknown_e;
        if( numberOfChanges == 0 ) {
            reaction->reactionType = MCGIDI_reactionType_scattering_e; }
        else {
            reaction->reactionType = MCGIDI_reactionType_nuclearIsomerTransmutation_e;
        }

/*
    Currently, these are not handled properly:
        MCGIDI_reactionType_nuclearLevelTransition_e
        MCGIDI_reactionType_atomic_e
*/
        break;
    }

    MCGIDI_reaction_setENDL_CSNumbers( smr, reaction );
    return( 0 );
}
/*
************************************************************
*/
static int MCGIDI_reaction_particleChanges( MCGIDI_POP *projectile, MCGIDI_POP *target, MCGIDI_productsInfo *productsInfo, int n1, int *particlesChanges ) {

    int projectileGlobalIndex = projectile->globalPoPsIndex, targetGlobalIndex = target->globalPoPsIndex, i1, i2 = 0;
    int gammaIndex = PoPs_particleIndex( "gamma" );

    if( projectileGlobalIndex != gammaIndex ) {
        for( i1 = 0; i1 < productsInfo->numberOfProducts; i1++ ) if( projectileGlobalIndex == productsInfo->productInfo[i1].globalPoPsIndex ) break;
        if( i1 == productsInfo->numberOfProducts ) particlesChanges[i2++] = projectileGlobalIndex;
    }

    for( i1 = 0; i1 < productsInfo->numberOfProducts; i1++ ) if( targetGlobalIndex == productsInfo->productInfo[i1].globalPoPsIndex ) break;
    if( i1 == productsInfo->numberOfProducts ) particlesChanges[i2++] = targetGlobalIndex;

    for( i1 = 0; i1 < productsInfo->numberOfProducts; i1++ ) {
        if( i2 == n1 ) break;
        if( /*(*/ projectileGlobalIndex == productsInfo->productInfo[i1].globalPoPsIndex /*)*/ ) continue;
        if( /*(*/ targetGlobalIndex == productsInfo->productInfo[i1].globalPoPsIndex /*)*/ ) continue;
        if( /*(*/ gammaIndex == productsInfo->productInfo[i1].globalPoPsIndex /*)*/ ) continue;
        particlesChanges[i2++] = productsInfo->productInfo[i1].globalPoPsIndex;
    }
    return( i2 );
}
/*
************************************************************
*/
static int MCGIDI_reaction_ParseDetermineReactionProducts( statusMessageReporting *smr, MCGIDI_POPs *pops, MCGIDI_outputChannel *outputChannel,
    MCGIDI_productsInfo *productsInfo, MCGIDI_reaction *reaction, double *finalQ, int level ) {
/*
*   This function determines all products that can be returned during sampling for this outputChannel. Note, products like 'U238_c' and
*   'U238_e3' are not returned during sampling as both are decay to the groud state (unless a meta-stable is encountered). 
*   Some examples for projectile 'n' and target 'U238' are:
*       outputChannel                       products returned during sampling.
*       'n + U238'                          n, U238
*       'n + U238 + gamma'                  n, U238, gamma
*       'n + U238_c'                        n, U238     (even if no gammas are give, return ground state of residual.
*       'n + (U238_e3 -> U238 + gamma)'     n, U238, gamma
*/
    int iProduct, nProducts = MCGIDI_outputChannel_numberOfProducts( outputChannel ), globalPoPsIndex, productIsTrackable;
    int twoBodyProductsWithData = 0;
    MCGIDI_product *product;
    MCGIDI_POP *residual;

    if( ( level == 0 ) && ( outputChannel->genre == MCGIDI_channelGenre_twoBody_e ) ) {
        for( iProduct = 0; iProduct < nProducts; iProduct++ ) {
            product = MCGIDI_outputChannel_getProductAtIndex( smr, outputChannel, iProduct );
            if( product->pop->globalPoPsIndex < 0 ) {
                twoBodyProductsWithData = -1; }
            else if( product->distribution.type == MCGIDI_distributionType_angular_e ) {
                if( twoBodyProductsWithData >= 0 ) twoBodyProductsWithData = 1;
            }
        }
    }
    if( twoBodyProductsWithData < 0 ) twoBodyProductsWithData = 0;
    *finalQ += MCGIDI_outputChannel_getQ_MeV( smr, outputChannel, 0 );
    for( iProduct = 0; iProduct < nProducts; iProduct++ ) {
        productIsTrackable = twoBodyProductsWithData;
        product = MCGIDI_outputChannel_getProductAtIndex( smr, outputChannel, iProduct );
        globalPoPsIndex = product->pop->globalPoPsIndex;
        if( ( product->distribution.type != MCGIDI_distributionType_none_e ) && ( product->distribution.type != MCGIDI_distributionType_unknown_e ) ) {
            productIsTrackable = 1;
            if( globalPoPsIndex < 0 ) {
                if( product->distribution.angular != NULL ) {
                    if( product->distribution.angular->type == MCGIDI_angularType_recoil ) productIsTrackable = 0;
                }
                if( productIsTrackable ) {
                    int len = (int) strlen( product->pop->name );

                    if( len > 2 ) {     /* Special case for continuum reactions with data for residual (e.g., n + U233 -> n + U233_c). */
                        if( ( product->pop->name[len-2] == '_' ) && ( product->pop->name[len-1] == 'c' ) ) {
                            for( residual = product->pop; residual->globalPoPsIndex < 0; residual = residual->parent ) ;
                            productIsTrackable = 1;
                            globalPoPsIndex = residual->globalPoPsIndex;
                        }
                    }
                    if( globalPoPsIndex < 0 ) {
                        smr_setReportError2( smr, smr_unknownID, 1, "product determination for '%s' cannot be determined", product->pop->name );
                        return( 1 );
                    }
                }
            }
        }
        if( productIsTrackable ) {
            if( MCGIDI_reaction_addReturnProduct( smr, productsInfo, globalPoPsIndex, product, reaction, 1 ) != 0 ) return( 1 ); }
        else {
            if( product->decayChannel.genre != MCGIDI_channelGenre_undefined_e ) {
                if( MCGIDI_reaction_ParseDetermineReactionProducts( smr, pops, &(product->decayChannel), productsInfo, reaction, finalQ, level + 1 ) != 0 ) return( 1 ); }
            else {
                *finalQ += product->pop->level_MeV;
                for( residual = product->pop; residual->globalPoPsIndex < 0; residual = residual->parent ) ;
                if( MCGIDI_reaction_addReturnProduct( smr, productsInfo, residual->globalPoPsIndex, product, reaction, 0 ) != 0 ) return( 1 );
                if( product->pop->numberOfGammaBranchs != 0 ) {
                    int gammaIndex = PoPs_particleIndex( "gamma" );
                    if( MCGIDI_reaction_addReturnProduct( smr, productsInfo, gammaIndex, NULL, reaction, 1 ) != 0 ) return( 1 );
                }
            }
        }
    }
    return( 0 );
}
/*
************************************************************
*/
static int MCGIDI_reaction_addReturnProduct( statusMessageReporting *smr, MCGIDI_productsInfo *productsInfo, int ID, MCGIDI_product *product, 
        MCGIDI_reaction *reaction, int transportable ) {

    int i1;
    enum MCGIDI_productMultiplicityType productMultiplicityType;

    MCGIDI_misc_updateTransportabilitiesMap2( reaction->transportabilities, ID, transportable );
    for( i1 = 0; i1 < productsInfo->numberOfProducts; i1++ ) {
        if( productsInfo->productInfo[i1].globalPoPsIndex == ID ) break;
    }
    if( i1 == productsInfo->numberOfProducts ) {
        if( productsInfo->numberOfProducts == productsInfo->numberOfAllocatedProducts ) {
            productsInfo->numberOfAllocatedProducts += 4;
            if( ( productsInfo->productInfo = (MCGIDI_productInfo *) smr_realloc2( smr, productsInfo->productInfo, 
                productsInfo->numberOfAllocatedProducts * sizeof( MCGIDI_productInfo ), "productsInfo->productInfo" ) ) == NULL ) return( 1 );
        }
        productsInfo->numberOfProducts++;
        productsInfo->productInfo[i1].globalPoPsIndex = ID;
        productsInfo->productInfo[i1].productMultiplicityType = MCGIDI_productMultiplicityType_unknown_e;
        productsInfo->productInfo[i1].multiplicity = 0;
        productsInfo->productInfo[i1].transportable = transportable;
    }
    if( product == NULL ) {
        productMultiplicityType = MCGIDI_productMultiplicityType_gammaBranching_e; }
    else {
        if( ( product->multiplicityVsEnergy != NULL ) || ( product->piecewiseMultiplicities != NULL ) ) {
            productMultiplicityType = MCGIDI_productMultiplicityType_energyDependent_e; }
        else {
            productsInfo->productInfo[i1].multiplicity += product->multiplicity;
            productMultiplicityType = MCGIDI_productMultiplicityType_integer_e;
        }
    }
    if( ( productsInfo->productInfo[i1].productMultiplicityType == MCGIDI_productMultiplicityType_unknown_e ) ||
        ( productsInfo->productInfo[i1].productMultiplicityType == productMultiplicityType ) ) {
        productsInfo->productInfo[i1].productMultiplicityType = productMultiplicityType; }
    else {
        productsInfo->productInfo[i1].productMultiplicityType = MCGIDI_productMultiplicityType_mixed_e;
    }
    return( 0 );
}
/*
************************************************************
*/
enum MCGIDI_reactionType MCGIDI_reaction_getReactionType( statusMessageReporting * /*smr*/, MCGIDI_reaction *reaction ) {

    return( reaction->reactionType );
}
/*
************************************************************
*/
MCGIDI_target_heated *MCGIDI_reaction_getTargetHeated( statusMessageReporting * /*smr*/, MCGIDI_reaction *reaction ) {

    return( reaction->target );
}
/*
************************************************************
*/
double MCGIDI_reaction_getProjectileMass_MeV( statusMessageReporting *smr, MCGIDI_reaction *reaction ) {

    return( MCGIDI_target_heated_getProjectileMass_MeV( smr, reaction->target ) );
}
/*
************************************************************
*/
double MCGIDI_reaction_getTargetMass_MeV( statusMessageReporting *smr, MCGIDI_reaction *reaction ) {

    return( MCGIDI_target_heated_getTargetMass_MeV( smr, reaction->target ) );
}
/*
************************************************************
*/
int MCGIDI_reaction_getDomain( statusMessageReporting * /*smr*/, MCGIDI_reaction *reaction, double *EMin, double *EMax ) {
/*
*   Return value
*       <  0    No cross section data.
*       == 0    Okay and EMin and EMax set.
*       >  0    error, EMin and EMax undefined.
*/

    if( !reaction->domainValuesPresent ) return( -1 );
    *EMin = reaction->EMin;
    *EMax = reaction->EMax;
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_reaction_fixDomains( statusMessageReporting * /*smr*/, MCGIDI_reaction *reaction, double EMin, double EMax, nfu_status *status ) {

    double lowerEps = 1e-14, upperEps = -1e-14;

    if( reaction->EMin == EMin ) lowerEps = 0.;
    if( reaction->EMax == EMax ) upperEps = 0.;
    if( ( lowerEps == 0. ) && ( upperEps == 0. ) ) return( 0 );

    *status = ptwXY_dullEdges( reaction->crossSection, lowerEps, upperEps, 1 );
    return( *status != nfu_Okay );
}
/*
************************************************************
*/
double MCGIDI_reaction_getCrossSectionAtE( statusMessageReporting *smr, MCGIDI_reaction *reaction, MCGIDI_quantitiesLookupModes &modes,
        bool sampling ) {

    double e_in = modes.getProjectileEnergy( ), xsec;

    if( modes.getCrossSectionMode( ) == MCGIDI_quantityLookupMode_pointwise ) {
        if( e_in < reaction->EMin ) e_in = reaction->EMin;
        if( e_in > reaction->EMax ) e_in = reaction->EMax;
        ptwXY_getValueAtX( reaction->crossSection, e_in, &xsec ); }
    else if( modes.getCrossSectionMode( ) == MCGIDI_quantityLookupMode_grouped ) {
        int index = modes.getGroupIndex( );
        double *xSecP = ptwX_getPointAtIndex( reaction->crossSectionGrouped, index );

        if( xSecP != NULL ) {
            xsec = *xSecP;
            if( sampling && ( index == reaction->thresholdGroupIndex ) ) xsec += reaction->thresholdGroupedDeltaCrossSection; }
        else {
            xsec = 0.;
            smr_setReportError2( smr, smr_unknownID, 1, "Invalid cross section group index %d", index );
        } }
    else {
        xsec = 0.;
    }
    return( xsec );
}
/*
************************************************************
*/
double MCGIDI_reaction_getFinalQ( statusMessageReporting * /*smr*/, MCGIDI_reaction *reaction, MCGIDI_quantitiesLookupModes &/*modes*/ ) {

    return( reaction->finalQ );
}
/*
************************************************************
*/
int MCGIDI_reaction_getENDF_MTNumber( MCGIDI_reaction *reaction ) {

    return( reaction->ENDF_MT );
}
/*
************************************************************
*/
int MCGIDI_reaction_getENDL_CSNumbers( MCGIDI_reaction *reaction, int *S ) {

    if( S != NULL ) *S = reaction->ENDL_S;
    return( reaction->ENDL_C );
}
/*
************************************************************
*/
static int MCGIDI_reaction_setENDL_CSNumbers( statusMessageReporting * /*smr*/, MCGIDI_reaction *reaction ) {

    int MT = MCGIDI_reaction_getENDF_MTNumber( reaction );
    int MT1_50ToC[] = { 1,   10,  -3,   -4,   -5,    0,    0,    0,    0,  -10,
                       32,    0,   0,    0,    0,   12,   13,   15,   15,   15,
                       15,   26,  36,   33,  -25,    0,  -27,   20,   27,  -30,
                        0,   22,  24,   25,  -35,  -36,   14,   15,    0,    0,
                       29,   16,   0,   17,   34,    0,    0,    0,    0 };
    int MT100_200ToC[] = { -101,   46,   40,   41,   42,   44,   45,   37, -109,    0,
                             18,   48, -113, -114,   19,   39,   47,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                              0, -152, -153, -154,   43, -156, -157,   23,   31, -160,
                           -161, -162, -163, -164, -165, -166, -167, -168, -169, -170,
                           -171, -172, -173, -174, -175, -176, -177, -178, -179, -180,
                           -181, -182, -183, -184, -185, -186, -187, -188,   28, -190,
                           -191, -192,   38, -194, -195, -196, -197, -198, -199, -200 };

    reaction->ENDL_C = 0;
    reaction->ENDL_S = 0;
    if( MT <= 0 ) return( 1 );
    if( MT > 891 ) return( 1 );
    if( MT < 50 ) {
        reaction->ENDL_C = MT1_50ToC[MT - 1]; }
    else if( MT <= 91 ) {
        reaction->ENDL_C = 11;
        if( MT != 91 ) reaction->ENDL_S = 1; }
    else if( ( MT > 100 ) && ( MT <= 200 ) ) {
        reaction->ENDL_C = MT100_200ToC[MT - 101]; }
    else if( ( MT == 452 ) || ( MT == 455 ) || ( MT == 456 ) || ( MT == 458 ) ) {
        reaction->ENDL_C = 15;
        if( MT == 455 ) reaction->ENDL_S = 7; }
    else if( MT >= 600 ) {
        if( MT < 650 ) {
            reaction->ENDL_C = 40;
            if( MT != 649 ) reaction->ENDL_S = 1; }
        else if( MT < 700 ) {
            reaction->ENDL_C = 41;
            if( MT != 699 ) reaction->ENDL_S = 1; }
        else if( MT < 750 ) {
            reaction->ENDL_C = 42;
            if( MT != 749 ) reaction->ENDL_S = 1; }
        else if( MT < 800 ) {
            reaction->ENDL_C = 44;
            if( MT != 799 ) reaction->ENDL_S = 1; }
        else if( MT < 850 ) {
            reaction->ENDL_C = 45;
            if( MT != 849 ) reaction->ENDL_S = 1; }
        else if( ( MT >= 875 ) && ( MT <= 891 ) ) {
            reaction->ENDL_C = 12;
            if( MT != 891 ) reaction->ENDL_S = 1;
        }
    }
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_productsInfo *MCGIDI_reaction_getProductsInfo( MCGIDI_reaction *reaction ) {

    return( &(reaction->productsInfo) );
}
/*
************************************************************
*/
int MCGIDI_reaction_recast( statusMessageReporting *smr, MCGIDI_reaction *reaction, GIDI_settings & /*settings*/, 
        GIDI_settings_particle const *projectileSettings, double temperature_MeV, ptwXPoints *totalGroupedCrossSection ) {

    if( totalGroupedCrossSection != NULL ) {
        nfu_status status_nf;
        GIDI_settings_group group( projectileSettings->getGroup( ) );

        if( reaction->crossSectionGrouped != NULL ) reaction->crossSectionGrouped = ptwX_free( reaction->crossSectionGrouped );
        if( ( reaction->crossSectionGrouped = projectileSettings->groupFunction( smr, reaction->crossSection, temperature_MeV, 0 ) ) == NULL ) return( 1 );
        if( ( status_nf = ptwX_add_ptwX( totalGroupedCrossSection, reaction->crossSectionGrouped ) ) != nfu_Okay ) return( 1 );

        reaction->thresholdGroupDomain = reaction->thresholdGroupedDeltaCrossSection = 0.;
        reaction->thresholdGroupIndex = group.getGroupIndexFromEnergy( reaction->EMin, false );
        if( reaction->thresholdGroupIndex > -1 ) {
            reaction->thresholdGroupDomain = group[reaction->thresholdGroupIndex+1] - reaction->EMin;
            if( reaction->thresholdGroupDomain > 0 ) {
                                                            /* factor 2 for linear reject in bin but above threshold. */
                reaction->thresholdGroupedDeltaCrossSection = *ptwX_getPointAtIndex( reaction->crossSectionGrouped, reaction->thresholdGroupIndex ) *
                        ( 2 * ( group[reaction->thresholdGroupIndex+1] - group[reaction->thresholdGroupIndex] ) / reaction->thresholdGroupDomain - 1 );
            }
        }
    }
    return( 0 );
}

/*
*********************** productsInfo ***********************
*/
/*
************************************************************
*/
int MCGIDI_productsInfo_getNumberOfUniqueProducts( MCGIDI_productsInfo *productsInfo ) {

    return( productsInfo->numberOfProducts );
}
/*
************************************************************
*/
int MCGIDI_productsInfo_getPoPsIndexAtIndex( MCGIDI_productsInfo *productsInfo, int index ) {

    if( ( index < 0 ) || ( index >= productsInfo->numberOfProducts ) ) return( -1 );
    return( productsInfo->productInfo[index].globalPoPsIndex );
}
/*
************************************************************
*/
enum MCGIDI_productMultiplicityType MCGIDI_productsInfo_getMultiplicityTypeAtIndex( MCGIDI_productsInfo *productsInfo, int index ) {

    if( ( index < 0 ) || ( index >= productsInfo->numberOfProducts ) ) return( MCGIDI_productMultiplicityType_invalid_e );
    return( productsInfo->productInfo[index].productMultiplicityType );
}
/*
************************************************************
*/
int MCGIDI_productsInfo_getIntegerMultiplicityAtIndex( MCGIDI_productsInfo *productsInfo, int index ) {

    if( ( index < 0 ) || ( index >= productsInfo->numberOfProducts ) ) return( -1 );
    return( productsInfo->productInfo[index].multiplicity );
}
/*
************************************************************
*/
int MCGIDI_productsInfo_getTransportableAtIndex( MCGIDI_productsInfo *productsInfo, int index ) {

    if( ( index < 0 ) || ( index >= productsInfo->numberOfProducts ) ) return( -1 );
    return( productsInfo->productInfo[index].transportable );
}

#if defined __cplusplus
}
#endif

