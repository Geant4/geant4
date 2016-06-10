/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <map>
#include <string.h>
#include <cmath>

#include "MCGIDI.h"
#include "MCGIDI_misc.h"
#include <xDataTOM_importXML_private.h>

#if defined __cplusplus
    extern "C" {
namespace GIDI {
using namespace GIDI;
#endif
static int _MCGIDI_target_releaseAndReturnOne( statusMessageReporting *smr, MCGIDI_target *target );
#if defined __cplusplus
    }
    }
#endif
/*
************************************************************
*/

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

MCGIDI_target *MCGIDI_target_new( statusMessageReporting *smr ) {

    MCGIDI_target *target;

    if( ( target = (MCGIDI_target *) smr_malloc2( smr, sizeof( MCGIDI_target ), 0, "target" ) ) == NULL ) return( NULL );
    if( MCGIDI_target_initialize( smr, target ) ) target = MCGIDI_target_free( smr, target );
    return( target );
}
/*
************************************************************
*/
int MCGIDI_target_initialize( statusMessageReporting * /*smr*/, MCGIDI_target *target ) {

    memset( target, 0, sizeof( MCGIDI_target ) );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_target *MCGIDI_target_newRead( statusMessageReporting *smr, const char *fileName ) {

    MCGIDI_target *target;

    if( ( target = MCGIDI_target_new( smr ) ) == NULL ) return( NULL );
    if( MCGIDI_target_read( smr, target, fileName ) != 0 ) smr_freeMemory( (void **) &target );
    return( target );
}
/*
************************************************************
*/
int MCGIDI_target_readFromMapViaPoPIDs( statusMessageReporting *smr, MCGIDI_target *target, MCGIDI_map *map, const char *evaluation, 
        int projectile_PoPID, int target_PoPID ) {

    char *targetPath;

    if( ( targetPath = MCGIDI_map_findTargetViaPoPIDs( smr, map, evaluation, projectile_PoPID, target_PoPID ) ) == NULL ) return( 1 );
    return( MCGIDI_target_read( smr, target, targetPath ) );
}
/*
************************************************************
*/
int MCGIDI_target_readFromMap( statusMessageReporting *smr, MCGIDI_target *target, MCGIDI_map *map, const char *evaluation, const char *projectileName, 
        const char *targetName ) {

    char *targetPath;

    if( ( targetPath = MCGIDI_map_findTarget( smr, map, evaluation, projectileName, targetName ) ) == NULL ) return( 1 );
    return( MCGIDI_target_read( smr, target, targetPath ) );
}
/*
************************************************************
*/
MCGIDI_target *MCGIDI_target_newReadFromMapViaPoPIDs( statusMessageReporting *smr, MCGIDI_map *map, const char *evaluation, 
        int projectile_PoPID, int target_PoPID ) {

    char *targetPath;
    MCGIDI_target *target;

    if( ( targetPath = MCGIDI_map_findTargetViaPoPIDs( smr, map, evaluation, projectile_PoPID, target_PoPID ) ) == NULL ) return( NULL );
    target = MCGIDI_target_newRead( smr, targetPath );
    smr_freeMemory( (void **) &targetPath );
    return( target );
}
/*
************************************************************
*/
MCGIDI_target *MCGIDI_target_newReadFromMap( statusMessageReporting *smr, MCGIDI_map *map, const char *evaluation, const char *projectileName, 
        const char *targetName ) {

    char *targetPath;
    MCGIDI_target *target;

    targetPath = MCGIDI_map_findTarget( smr, map, evaluation, projectileName, targetName );
    if( targetPath == NULL ) return( NULL );
    target = MCGIDI_target_newRead( smr, targetPath );
    smr_freeMemory( (void **) &targetPath );
    return( target );
}
/*
************************************************************
*/
MCGIDI_target *MCGIDI_target_free( statusMessageReporting *smr, MCGIDI_target *target ) {

    MCGIDI_target_release( smr, target );
    smr_freeMemory( (void **) &target );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_target_release( statusMessageReporting *smr, MCGIDI_target *target ) {

    int i;

    smr_freeMemory( (void **) &(target->path) );
    smr_freeMemory( (void **) &(target->absPath) );
    xDataTOMAL_release( &(target->attributes) );
    for( i = 0; i < target->nHeatedTargets; i++ ) {
        smr_freeMemory( (void **) &(target->heatedTargets[i].path) );
        smr_freeMemory( (void **) &(target->heatedTargets[i].contents) );
        if( target->heatedTargets[i].heatedTarget != NULL ) MCGIDI_target_heated_free( smr, target->heatedTargets[i].heatedTarget );
    }
    smr_freeMemory( (void **) &(target->heatedTargets) );
    smr_freeMemory( (void **) &(target->readHeatedTargets) );
    MCGIDI_target_initialize( smr, target );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_target_read( statusMessageReporting *smr, MCGIDI_target *target, const char *fileName ) {
/*
*   If a target has already been read into this target, user must have called MCGIDI_target_release before calling this routine.
*   Otherwise, there will be memory leaks.
*/
    xDataXML_document *doc;
    xDataXML_element *element, *child;
    int i, iHeated, nHeated = 0, status = 1;
    double temperature;
    /* char *pReturnValue; */
    char const *version, *contents;

    MCGIDI_target_initialize( smr, target );
    if( ( target->path = smr_allocateCopyString2( smr, fileName, "path" ) ) == NULL ) return( status );
    if( ( target->absPath = MCGIDI_misc_getAbsPath( smr, fileName ) ) == NULL ) return( _MCGIDI_target_releaseAndReturnOne( smr, target ) );
    if( ( doc = xDataXML_importFile2( smr, fileName ) ) == NULL ) return( _MCGIDI_target_releaseAndReturnOne( smr, target ) );
    element = xDataXML_getDocumentsElement( doc );
    if( strcmp( element->name, "xTarget" ) != 0 ) {
        MCGIDI_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "input file's top element must be xTarget and not %s", element->name ); }
    else {
        status = 0;
        if( ( version = xDataXML_getAttributesValueInElement( element, "version" ) ) == NULL ) {
            smr_setReportError2( smr, smr_unknownID, 1, "version attribute missing from element '%s'", element->name );
            status = 1; }
        else {
            if( strcmp( version, "xMCProcess 0.1" ) != 0 ) {
                smr_setReportError2( smr, smr_unknownID, 1, "Unsupported version '%s' for element %s", version, element->name );
                status = 1;
            }
        }
        if( status == 0 ) {
          /*  pReturnValue = ( MCGIDI_misc_copyXMLAttributesToTOM( smr, &(target->attributes), &(element->attributes) ) ) ? NULL : target->path; */
            for( nHeated = 0, child = xDataXML_getFirstElement( element ); child != NULL; nHeated++, child = xDataXML_getNextElement( child ) ) {
                if( strcmp( child->name, "target" ) != 0 ) {
                    MCGIDI_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "element can only have target sub-elements%s", 
                                element->name );
                    status = 1;
                    break;
                }
            }
        }
        if( status == 0 ) {
            if( ( target->heatedTargets = (MCGIDI_target_heated_info *) smr_malloc2( smr, nHeated * sizeof( MCGIDI_target_heated_info ), 1, "heatedTargets" ) ) == NULL ) {
                status = 1; }
            else {
                if( ( target->readHeatedTargets = (MCGIDI_target_heated_info **) smr_malloc2( smr, nHeated * sizeof( MCGIDI_target_heated_info * ), 1, "heatedTargets" ) ) == NULL ) 
                    status = 1;
            }
            for( nHeated = 0, child = xDataXML_getFirstElement( element ); ( status == 0 ) && ( child != NULL ); 
                    nHeated++, child = xDataXML_getNextElement( child ) ) {
                if( ( i = xDataXML_convertAttributeToDouble( smr, child, "temperature", &temperature, 1 ) ) != 0 ) {
                    if( i > 0 ) smr_setReportError2p( smr, smr_unknownID, 1, "target does not have a temperature attribute" );
                    status = 1;
                    break;
                }
                for( iHeated = 0; iHeated < nHeated; iHeated++ ) if( target->heatedTargets[iHeated].temperature > temperature ) break;
                if( iHeated < nHeated ) for( i = nHeated; i >= iHeated; i-- ) target->heatedTargets[i+1] = target->heatedTargets[i];
                target->heatedTargets[iHeated].temperature = temperature;
                target->heatedTargets[iHeated].path = NULL;
                target->heatedTargets[iHeated].contents = NULL;
                target->heatedTargets[iHeated].heatedTarget = NULL;
                if( ( contents = xDataXML_getAttributesValueInElement( child, "contents" ) ) != NULL ) {
                    if( ( target->heatedTargets[iHeated].contents = smr_allocateCopyString2( smr, contents, "contents" ) ) == NULL ) {
                        status = 1;
                        break;
                    }
                }
                if( ( contents = xDataXML_getAttributesValueInElement( child, "file" ) ) == NULL ) {
                    status = 1;
                    break;
                }
                if( ( target->heatedTargets[iHeated].path = (char *) smr_malloc2(smr, strlen( target->absPath ) + strlen( contents ) + 2, 0, "path") ) == NULL) {
                    status = 1;
                    break;
                }
                strcpy( target->heatedTargets[iHeated].path, target->absPath );
                *strrchr( target->heatedTargets[iHeated].path, '/' ) = 0;
                strcat( target->heatedTargets[iHeated].path, "/" );
                strcat( target->heatedTargets[iHeated].path, contents );
                target->nHeatedTargets++;
            }
        }
    }
    xDataXML_freeDoc( smr, doc );
    if( status == 0 ) {
        for( i = 0; i < nHeated; i++ ) target->heatedTargets[i].ordinal = i;
        for( i = 0; i < nHeated; i++ ) if( target->heatedTargets[i].contents == NULL ) break;
        if( i == nHeated ) i = 0;                                           /* All heated targets are crossSection only. */
        if( MCGIDI_target_readHeatedTarget( smr, target, i ) == 0 ) {
            target->baseHeatedTarget = target->heatedTargets[i].heatedTarget; }
        else {
            MCGIDI_target_release( NULL, target );
            status = 1;
        } }
    else {
        MCGIDI_target_release( smr, target );
    }
    return( status );
}
/*
************************************************************
*/
char const *MCGIDI_target_getAttributesValue( statusMessageReporting * /*smr*/, MCGIDI_target *target, char const *name ) {

    return( xDataTOMAL_getAttributesValue( &(target->attributes), name ) );
}
/*
************************************************************
*/
int MCGIDI_target_getTemperatures( statusMessageReporting * /*smr*/, MCGIDI_target *target, double *temperatures ) {

    int i;

    if( temperatures != NULL ) for( i = 0; i < target->nHeatedTargets; i++ ) temperatures[i] = target->heatedTargets[i].temperature;
    return( target->nHeatedTargets );
}
/*
************************************************************
*/
int MCGIDI_target_readHeatedTarget( statusMessageReporting *smr, MCGIDI_target *target, int index ) {

    int i;

    if( ( index < 0 ) || ( index >= target->nHeatedTargets ) ) {
        smr_setReportError2( smr, smr_unknownID, 1, "temperature index = %d out of range (0 <= index < %d", index, target->nHeatedTargets );
        return( -1 );
    }
    if( target->heatedTargets[index].heatedTarget != NULL ) return( 1 );
    if( ( target->heatedTargets[index].heatedTarget = MCGIDI_target_heated_newRead( smr, target->heatedTargets[index].path ) ) != NULL ) {
        target->projectilePOP = target->heatedTargets[index].heatedTarget->projectilePOP;
        target->targetPOP = target->heatedTargets[index].heatedTarget->targetPOP;
        if( target->heatedTargets[index].heatedTarget != NULL ) {
            target->heatedTargets[index].heatedTarget->ordinal = target->heatedTargets[index].ordinal;
            for( i = target->nReadHeatedTargets; i > 0; i-- ) {
                if( target->readHeatedTargets[i-1]->temperature < target->heatedTargets[index].temperature ) break;
                target->readHeatedTargets[i] = target->readHeatedTargets[i-1];
            }
            target->readHeatedTargets[i] = &(target->heatedTargets[i]);
            target->nReadHeatedTargets++;
        }
    }
    return( ( target->heatedTargets[index].heatedTarget == NULL ? -1 : 0 ) );
}
/*
************************************************************
*/
MCGIDI_target_heated *MCGIDI_target_getHeatedTargetAtIndex_ReadIfNeeded( statusMessageReporting *smr, MCGIDI_target *target, int index ) {

    if( ( index < 0 ) || ( index >= target->nHeatedTargets ) ) {
        smr_setReportError2( smr, smr_unknownID, 1, "temperature index = %d out of range (0 <= index < %d", index, target->nHeatedTargets );
        return( NULL );
    }
    if( target->heatedTargets[index].heatedTarget == NULL ) MCGIDI_target_readHeatedTarget( smr, target, index );
    return( target->heatedTargets[index].heatedTarget );
}
/*
************************************************************
*/
MCGIDI_target_heated *MCGIDI_target_getHeatedTargetAtTIndex( statusMessageReporting *smr, MCGIDI_target *target, int index ) {

    if( ( index < 0 ) || ( index >= target->nHeatedTargets ) ) {
        smr_setReportError2( smr, smr_unknownID, 1, "temperature index = %d out of range (0 <= index < %d", index, target->nHeatedTargets );
        return( NULL );
    }
    if( target->heatedTargets[index].heatedTarget == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "temperature index = %d not read in", index );
        return( NULL );
    }
    return( target->heatedTargets[index].heatedTarget );
}
/*
************************************************************
*/
int MCGIDI_target_numberOfReactions( statusMessageReporting *smr, MCGIDI_target *target ) {

    return( MCGIDI_target_heated_numberOfReactions( smr, target->baseHeatedTarget ) );
}
/*
************************************************************
*/
enum MCGIDI_reactionType MCGIDI_target_getReactionTypeAtIndex( statusMessageReporting *smr, MCGIDI_target *target, int index ) {

    MCGIDI_reaction *reaction = MCGIDI_target_heated_getReactionAtIndex_smr( smr, target->baseHeatedTarget, index );

    if( reaction == NULL ) return( MCGIDI_reactionType_unknown_e );
    return( MCGIDI_reaction_getReactionType( smr, reaction ) );
}
/*
************************************************************
*/
MCGIDI_reaction *MCGIDI_target_getReactionAtIndex( MCGIDI_target *target, int index ) {

    return( MCGIDI_target_heated_getReactionAtIndex( target->baseHeatedTarget, index ) );
}
/*
************************************************************
*/
MCGIDI_reaction *MCGIDI_target_getReactionAtIndex_smr( statusMessageReporting *smr, MCGIDI_target *target, int index ) {

    return( MCGIDI_target_heated_getReactionAtIndex_smr( smr, target->baseHeatedTarget, index ) );
}
/*
************************************************************
*/
int MCGIDI_target_numberOfProductionReactions( statusMessageReporting * /*smr*/, MCGIDI_target * /*target*/ ) {

#if 0
    return( MCGIDI_target_heated_numberOfProductionReactions( smr, target->baseHeatedTarget ) );
#endif
return( 0 );
}

/*
************************************************************
*/
double MCGIDI_target_getTotalCrossSectionAtTAndE( statusMessageReporting *smr, MCGIDI_target *target, MCGIDI_quantitiesLookupModes &modes,
        bool sampling ) {

    int i;
    double xsec = 0., xsec1, xsec2, temperature = modes.getTemperature( );

    for( i = 0; i < target->nReadHeatedTargets; i++ ) if( target->readHeatedTargets[i]->temperature > temperature ) break;
    if( i == 0 ) {
        xsec = MCGIDI_target_heated_getTotalCrossSectionAtE( smr, target->readHeatedTargets[0]->heatedTarget, modes, sampling ); }
    else if( i == target->nReadHeatedTargets ) {
        xsec = MCGIDI_target_heated_getTotalCrossSectionAtE( smr, target->readHeatedTargets[i-1]->heatedTarget, modes, sampling ); }
    else {
        xsec1 = MCGIDI_target_heated_getTotalCrossSectionAtE( smr, target->readHeatedTargets[i-1]->heatedTarget, modes, sampling );
        xsec2 = MCGIDI_target_heated_getTotalCrossSectionAtE( smr, target->readHeatedTargets[i  ]->heatedTarget, modes, sampling );
        xsec = ( ( target->readHeatedTargets[i]->temperature - temperature ) * xsec1 + 
                 ( temperature - target->readHeatedTargets[i-1]->temperature ) * xsec2 ) / 
                 ( target->readHeatedTargets[i]->temperature - target->readHeatedTargets[i-1]->temperature );
    }

    return( xsec );
}
/*
************************************************************
*/
int MCGIDI_target_getDomain( statusMessageReporting *smr, MCGIDI_target *target, double *EMin, double *EMax ) {

    int ir, nr = MCGIDI_target_numberOfReactions( smr, target );
    double EMin_, EMax_;

    for( ir = 0; ir < nr; ir++ ) {
        MCGIDI_target_heated_getReactionsDomain( smr, target->baseHeatedTarget, ir, &EMin_, &EMax_ );
        if( ir == 0 ) {
            *EMin = EMin_;
            *EMax = EMax_; }
        else {
            if( *EMin > EMin_ ) *EMin = EMin_;
            if( *EMax < EMax_ ) *EMax = EMax_;
        }
    }
    return( 0 );
}
/*
************************************************************
*/
double MCGIDI_target_getIndexReactionCrossSectionAtE( statusMessageReporting *smr, MCGIDI_target *target, int index, MCGIDI_quantitiesLookupModes &modes,
        bool sampling ) {

    int i;
    double xsec = 0., xsec1, xsec2, temperature = modes.getTemperature( );

    for( i = 0; i < target->nReadHeatedTargets; i++ ) if( target->readHeatedTargets[i]->temperature > temperature ) break;
    if( i == 0 ) {
        xsec = MCGIDI_target_heated_getIndexReactionCrossSectionAtE( smr, target->readHeatedTargets[0]->heatedTarget, index, modes, sampling ); }
    else if( i == target->nReadHeatedTargets ) {
        xsec = MCGIDI_target_heated_getIndexReactionCrossSectionAtE( smr, target->readHeatedTargets[i-1]->heatedTarget, index, modes, sampling ); }
    else {
        xsec1 = MCGIDI_target_heated_getIndexReactionCrossSectionAtE(smr, target->readHeatedTargets[i-1]->heatedTarget, index, modes, sampling );
        xsec2 = MCGIDI_target_heated_getIndexReactionCrossSectionAtE(smr, target->readHeatedTargets[i  ]->heatedTarget, index, modes, sampling );
        xsec = ( ( target->readHeatedTargets[i]->temperature - temperature ) * xsec1 + 
                 ( temperature - target->readHeatedTargets[i-1]->temperature ) * xsec2 ) / 
                 ( target->readHeatedTargets[i]->temperature - target->readHeatedTargets[i-1]->temperature );
    }

    return( xsec );
}
/*
************************************************************
*/
int MCGIDI_target_sampleReaction( statusMessageReporting *smr, MCGIDI_target *target, MCGIDI_quantitiesLookupModes &modes, double totalXSec, 
        double (*userrng)( void * ), void *rngState ) {

    int ir, nr = MCGIDI_target_numberOfReactions( smr, target );
    double rngValue = (*userrng)( rngState );
    double cumm_xsec = 0., r_xsec = rngValue * totalXSec;

    for( ir = 0; ir < nr; ir++ ) {
        cumm_xsec += MCGIDI_target_getIndexReactionCrossSectionAtE( smr, target, ir, modes, true );
        if( cumm_xsec >= r_xsec ) break;
    }
    if( ir == nr ) {
        if( ( totalXSec - cumm_xsec ) >= 1e-12 * totalXSec ) {
            smr_setReportError2( smr, smr_unknownID, 1, 
                "Failed to sample a reaction for temperature = %.12e, energy = %.12e, totalXSec = %16.e, rngValue = %16.e, r_xsec = %16.e, cumm_xsec = %16.e", 
                modes.getTemperature( ), modes.getProjectileEnergy( ), totalXSec, rngValue, r_xsec, cumm_xsec );
            return( -1 );
        }
        ir--;                               /* May not be correct but close. */
    }
    if( modes.getCrossSectionMode( ) == MCGIDI_quantityLookupMode_grouped ) {
        MCGIDI_reaction *reaction = MCGIDI_target_heated_getReactionAtIndex( target->baseHeatedTarget, ir );

        if( modes.getGroupIndex( ) == reaction->thresholdGroupIndex ) {
            double dEnergy = modes.getProjectileEnergy( ) - reaction->EMin;

            if( dEnergy  <= 0 ) return( MCGIDI_nullReaction );
            if( ( (*userrng)( rngState ) * reaction->thresholdGroupDomain ) > dEnergy ) return( MCGIDI_nullReaction );
        }
    }
    return( ir );
}
/*
************************************************************
*/
int MCGIDI_target_sampleNullReactionProductsAtE( statusMessageReporting *smr, MCGIDI_target *target,
    MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo, MCGIDI_sampledProductsDatas *productDatas ) {

    MCGIDI_sampledProductsData productData;

    productData.isVelocity = decaySamplingInfo->isVelocity;
    productData.pop = target->projectilePOP;
    productData.kineticEnergy = modes.getProjectileEnergy( );
    productData.px_vx = 0.;
    productData.py_vy = 0.;
    productData.pz_vz = std::sqrt( productData.kineticEnergy * ( productData.kineticEnergy + 2. * productData.pop->mass_MeV ) );
    if( productData.isVelocity ) productData.pz_vz *=
            MCGIDI_speedOfLight_cm_sec / std::sqrt( productData.pz_vz * productData.pz_vz + productData.pop->mass_MeV * productData.pop->mass_MeV );
    productData.delayedNeutronIndex = 0;
    productData.delayedNeutronRate = 0.;
    productData.birthTimeSec = 0;

    productDatas->numberOfProducts = 0;
    MCGIDI_sampledProducts_addProduct( smr, productDatas, &productData );
    return( productDatas->numberOfProducts );
}
/*
************************************************************
*/
int MCGIDI_target_sampleIndexReactionProductsAtE( statusMessageReporting *smr, MCGIDI_target *target, int index,
        MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo, MCGIDI_sampledProductsDatas *productData ) {

    return( MCGIDI_target_heated_sampleIndexReactionProductsAtE( smr, target->baseHeatedTarget, index, modes, decaySamplingInfo, productData ) );
}
/*
************************************************************
*/
double MCGIDI_target_getIndexReactionFinalQ( statusMessageReporting *smr, MCGIDI_target *target, int index, 
        MCGIDI_quantitiesLookupModes &modes ) {

    return( MCGIDI_target_heated_getIndexReactionFinalQ( smr, target->baseHeatedTarget, index, modes ) );
}
/*
************************************************************
*/
std::map<int, enum MCGIDI_transportability> const *MCGIDI_target_getUniqueProducts( statusMessageReporting *smr, MCGIDI_target *target ) {

    return( MCGIDI_target_heated_getUniqueProducts( smr, target->baseHeatedTarget ) );
}
/*
************************************************************
*/
int MCGIDI_target_recast( statusMessageReporting *smr, MCGIDI_target *target, GIDI_settings &settings ) {

    int i1, status = 0;

    for( i1 = 0; i1 < target->nReadHeatedTargets; i1++ ) {
        if( ( status = MCGIDI_target_heated_recast( smr, target->readHeatedTargets[i1]->heatedTarget, settings ) ) != 0 ) break;
    }
    return( status );
}
/*
************************************************************
*/
static int _MCGIDI_target_releaseAndReturnOne( statusMessageReporting *smr, MCGIDI_target *target ) {

    MCGIDI_target_release( smr, target );
    return( 1 );
}

#if defined __cplusplus
}
#endif
