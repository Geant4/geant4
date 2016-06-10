/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <map>
#include <string.h>
#include <cmath>

#include <xDataTOM.h>
#include "MCGIDI.h"
#include "MCGIDI_misc.h"
#include "MCGIDI_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static int MCGIDI_target_heated_parsePOPs( statusMessageReporting *smr, MCGIDI_target_heated *target, xDataTOM_element *element,
    xDataTOM_element *particleAliases );
static int MCGIDI_target_heated_parseParticle( statusMessageReporting *smr, MCGIDI_target_heated *target, xDataTOM_element *element,
    xDataTOM_element *particleAliases );
static int MCGIDI_target_heated_parseParticleLevel( statusMessageReporting *smr, MCGIDI_target_heated *target, xDataTOM_element *element, MCGIDI_POP *parent,
    double mass_MeV, xDataTOM_element *particleAliases );
static int MCGIDI_target_heated_parseParticleGammas( statusMessageReporting *smr, MCGIDI_target_heated *target, xDataTOM_element *element, char const *name );
static int MCGIDI_target_heated_parseReaction( statusMessageReporting *smr, xDataTOM_element *child, MCGIDI_target_heated *target, 
    MCGIDI_POPs *pops, MCGIDI_reaction *reaction );
/*
************************************************************
*/
MCGIDI_target_heated *MCGIDI_target_heated_new( statusMessageReporting *smr ) {

    MCGIDI_target_heated *target;

    if( ( target = (MCGIDI_target_heated *) smr_malloc2( smr, sizeof( MCGIDI_target_heated ), 0, "target" ) ) == NULL ) return( NULL );
    if( MCGIDI_target_heated_initialize( smr, target ) ) target = (MCGIDI_target_heated *) smr_freeMemory( (void **) &target );
    return( target );
}
/*
************************************************************
*/
int MCGIDI_target_heated_initialize( statusMessageReporting *smr, MCGIDI_target_heated *target ) {

    memset( target, 0, sizeof( MCGIDI_target_heated ) );
    MCGIDI_POPs_initial( smr, &(target->pops), 100 );
    target->transportabilities = new transportabilitiesMap( );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_target_heated *MCGIDI_target_heated_newRead( statusMessageReporting *smr, const char *fileName ) {

    MCGIDI_target_heated *target;

    if( ( target = MCGIDI_target_heated_new( smr ) ) == NULL ) return( NULL );
    if( MCGIDI_target_heated_read( smr, target, fileName ) != 0 ) target = (MCGIDI_target_heated *) smr_freeMemory( (void **) &target );
    return( target );
}
/*
************************************************************
*/
MCGIDI_target_heated *MCGIDI_target_heated_free( statusMessageReporting *smr, MCGIDI_target_heated *target ) {

    MCGIDI_target_heated_release( smr, target );
    smr_freeMemory( (void **) &target );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_target_heated_release( statusMessageReporting * /*smr*/, MCGIDI_target_heated *target ) {

    int ir;

    ptwXY_free( target->crossSection );
    ptwX_free( target->crossSectionGrouped );
    ptwX_free( target->crossSectionGroupedForSampling );
    for( ir = 0; ir < target->numberOfReactions; ir++ ) MCGIDI_reaction_release( NULL, &(target->reactions[ir]) );
    smr_freeMemory( (void **) &(target->reactions) );
    MCGIDI_POPs_release( &(target->pops) );
    smr_freeMemory( (void **) &(target->path) );
    smr_freeMemory( (void **) &(target->absPath) );
    xDataTOMAL_release( &(target->attributes) );
    delete target->transportabilities;
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_target_heated_read( statusMessageReporting *smr, MCGIDI_target_heated *target, const char *fileName ) {
/*
*   If a target has already been read into this target, user must have called MCGIDI_target_heated_release before calling this routine.
*   Otherwise, there will be memory leaks.
*/
    int n, ir;
    xDataTOM_TOM *doc = NULL;
    xDataTOM_element *element, *child, *particles, *particleAliases;
    char const *name, *version, *temperatureStr;
    char *e1;
    MCGIDI_reaction *reaction;
    double crossSectionInit[4] = { 0., 0., 0., 0., };
    nfu_status status;
    ptwXYPoints *crossSection;
int subtag1_Notice = 0;

    if( ( target->path = smr_allocateCopyString2( smr, fileName, "path" ) ) == NULL ) goto err;
    if( ( target->absPath = xDataTOMMisc_getAbsPath( smr, fileName ) ) == NULL ) goto err;
    if( ( doc = xDataTOM_importFile( smr, fileName ) ) == NULL ) goto err;
    element = xDataTOM_getDocumentsElement( doc );
    if( ( version = xDataTOM_getAttributesValueInElement( element, "version" ) ) == NULL ) {
            smr_setReportError2( smr, smr_unknownID, 1, "version attribute missing from element '%s'", element->name );
            goto err; }
    else {
        if( strcmp( version, "GND 1.3" ) != 0 ) {
            smr_setReportError2( smr, smr_unknownID, 1, "Unsupported version '%s' for element %s", version, element->name );
            goto err;
        }
    }
    if( strcmp( element->name, "reactionSuite" ) != 0 ) {
        smr_setReportError2( smr, smr_unknownID, 1, "input file's top element must be reactionSuite and not %s", element->name );
        goto err; }
    else {
        xDataTOMAL_copyAttributionList( smr, &(target->attributes), &(element->attributes) );
        particleAliases = xDataTOME_getOneElementByName( smr, element, "aliases", 0 );
        if( ( particles = xDataTOME_getOneElementByName( smr, element, "particles", 1 ) ) == NULL ) goto err;
        if( MCGIDI_target_heated_parsePOPs( smr, target, particles, particleAliases ) != 0 ) goto err;

        if( ( temperatureStr = MCGIDI_misc_pointerToTOMAttributeIfAllOk3( smr, target->absPath, 1, &(target->attributes), "temperature" ) ) == NULL ) goto err;
        target->temperature_MeV = strtod( temperatureStr, &e1 );
        while( isspace( *e1 ) ) ++e1; // Loop checking, 11.06.2015, T. Koi
        target->temperature_MeV *= MCGIDI_misc_getUnitConversionFactor( smr, e1, "MeV/k" );
        if( !smr_isOk( smr ) ) goto err;

        if( ( name = MCGIDI_misc_pointerToTOMAttributeIfAllOk3( smr, target->absPath, 1, &(target->attributes), "projectile" ) ) != NULL )
            target->projectilePOP = MCGIDI_POPs_findParticle( &(target->pops), name );
        if( !smr_isOk( smr ) ) goto err;

        if( ( name = MCGIDI_misc_pointerToTOMAttributeIfAllOk3( smr, target->absPath, 1, &(target->attributes), "target" ) ) != NULL )
        if( !smr_isOk( smr ) ) goto err;
            target->targetPOP = MCGIDI_POPs_findParticle( &(target->pops), name );

        n = xDataTOM_numberOfElementsByName( smr, element, "reaction" );
        if( n == 0 ) {
            smr_setReportError2( smr, smr_unknownID, 1, "target does not have any reactions: file = '%s'", fileName );
            goto err;
        }
        if( ( target->reactions = (MCGIDI_reaction *) smr_malloc2( smr, n * sizeof( MCGIDI_reaction ), 1, "target->reactions" ) ) == NULL ) goto err;

        for( ir = 0, child = xDataTOME_getFirstElement( element ); child != NULL; child = xDataTOME_getNextElement( child ) ) {
            if( strcmp( child->name, "particles" ) == 0 ) continue;
            if( strcmp( child->name, "styles" ) == 0 ) continue;
            if( strcmp( child->name, "documentations" ) == 0 ) continue;
            if( strcmp( child->name, "resonances" ) == 0 ) continue;
            if( strcmp( child->name, "summedReaction" ) == 0 ) continue;
            if( strcmp( child->name, "fissionComponent" ) == 0 ) continue;
            if( strcmp( child->name, "reaction" ) == 0 ) {
                double EMin, EMax;

                reaction =  &(target->reactions[ir]);
                if( MCGIDI_target_heated_parseReaction( smr, child, target, &(target->pops), reaction ) ) goto err;
                if( MCGIDI_reaction_getDomain( smr, reaction, &EMin, &EMax ) ) goto err;
                if( ir == 0 ) { target->EMin = EMin; target->EMax = EMax; }
                if( EMin < target->EMin ) target->EMin = EMin;
                if( EMax > target->EMax ) target->EMax = EMax;
                for( transportabilitiesMap::const_iterator iter = reaction->transportabilities->begin( ); 
                            iter != reaction->transportabilities->end( ); ++iter ) {
                    MCGIDI_misc_updateTransportabilitiesMap( target->transportabilities, iter->first, iter->second );
                }
                ir++; }
            else if( strcmp( child->name, "production" ) == 0 ) {
                continue; }
            else if( strcmp( child->name, "aliases" ) == 0 ) {
                continue; }
            else if( strcmp( child->name, "partialGammaProduction" ) == 0 ) {
                if( subtag1_Notice == 0 ) printf( "Unsupported reactionSuite sub-tag = '%s'\n", child->name );
                subtag1_Notice++; }
            else {
                printf( "Unsupported reactionSuite sub-tag = '%s'\n", child->name );
            }
        }
        crossSectionInit[0] = target->EMin;
        crossSectionInit[2] = target->EMax;
        if( ( target->crossSection = ptwXY_create( ptwXY_interpolationLinLin, NULL, 2., 1e-3, 2, 10, 2, crossSectionInit, &status, 0 ) ) == NULL ) {
            smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_create err = %d: %s\n", status, nfu_statusMessage( status ) );
            goto err;
        }
        for( ir = 0; ir < target->numberOfReactions; ir++ ) {
            reaction =  &(target->reactions[ir]);
            if( MCGIDI_reaction_fixDomains( smr, reaction, target->EMin, target->EMax, &status ) ) {
                smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_add_ptwXY err = %d: %s\n", status, nfu_statusMessage( status ) );
                goto err;
            }
            if( ( crossSection = ptwXY_add_ptwXY( target->crossSection, reaction->crossSection, &status ) ) == NULL ) {
                smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_add_ptwXY err = %d: %s\n", status, nfu_statusMessage( status ) );
                goto err;
            }
            target->crossSection = ptwXY_free( target->crossSection );
            target->crossSection = crossSection;
        }
    }
    xDataTOM_freeTOM( smr, &doc );
    return( 0 );

err:
    smr_setReportError2( smr, smr_unknownID, 1, "Sub-error while reading file '%s'", fileName );
    if( doc != NULL ) xDataTOM_freeTOM( smr, &doc );
    MCGIDI_target_heated_release( smr, target );
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_target_heated_parsePOPs( statusMessageReporting *smr, MCGIDI_target_heated *target, xDataTOM_element *element,
        xDataTOM_element *particleAliases ) {

    xDataTOM_element *child;

    for( child = xDataTOME_getFirstElement( element ); child != NULL; child = xDataTOME_getNextElement( child ) ) {
        if( strcmp( child->name, "particle" ) ) {
            smr_setReportError2( smr, smr_unknownID, 1, "invalid element '%s' in %s", child->name, element->name );
            goto err;
        }
        if( MCGIDI_target_heated_parseParticle( smr, target, child, particleAliases ) ) goto err;
    }
    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_target_heated_parseParticle( statusMessageReporting *smr, MCGIDI_target_heated *target, xDataTOM_element *element,
        xDataTOM_element *particleAliases ) {
/*
    This routine, MCGIDI_target_heated_parseParticleLevel, MCGIDI_target_heated_parseParticleGammas handle the parsing of a
    particle which can have one of the following three forms:

1)
    <particle name="Sc47" genre="nucleus" mass="46.9524440292728 amu"/>
2)
    <particle name="Sc46" genre="nucleus" mass="45.9551770270807 amu">
        <level name="Sc46_e0" label="0" energy="0 eV"/>
        <level name="Sc46_e1" label="1" energy="142500 eV"/></particle>
3)
    <particle name="S36" genre="nucleus" mass="35.9669735654569 amu">
        <level name="S36_e0" label="0" energy="0 eV" spin="0"/>
        <level name="S36_e1" label="1" energy="3.291e6 eV">
            <gamma finalLevel="S36_e0" probability="1.0"/></level>
        <level name="S36_e2" label="2" energy="3.346e6 eV">
            <gamma finalLevel="S36_e0" probability="1.0"/></level>
        <level name="S36_e3" label="3" energy="4191990 eV">
            <gamma finalLevel="S36_e1" probability="1.0"/></level>
        <level name="S36_e4" label="4" energy="4522990 eV">
            <gamma finalLevel="S36_e1" probability="0.248"/>
            <gamma finalLevel="S36_e0" probability="0.752"/></level>
        <level name="S36_e5" label="5" energy="4574990 eV">
            <gamma finalLevel="S36_e1" probability="1.0"/></level>
        <level name="S36_c" label="c" energy="u:4999990 eV"/></particle>
*/
    int globalParticle = 1;
    char const *name = NULL, *mass = NULL;      /* Do not free name or mass, do not own them. */
    double mass_MeV;
    xDataTOM_element *child;
    MCGIDI_POP *pop;

    if( ( name = xDataTOM_getAttributesValueInElement( element, "name" ) ) == NULL ) {
        smr_setReportError2p( smr, smr_unknownID, 1, "particle missing name attribute" );
        goto err;
    }
    if( ( mass = xDataTOM_getAttributesValueInElement( element, "mass" ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "particle '%s' missing mass attribute", name );
        goto err;
    }
    if( MCGIDI_misc_PQUStringToDouble( smr, mass, "amu", MCGIDI_AMU2MeV, &mass_MeV ) ) goto err;
    if( ( pop = MCGIDI_POPs_addParticleIfNeeded( smr, &(target->pops), name, mass_MeV, 0., NULL, globalParticle ) ) == NULL ) goto err;

    for( child = xDataTOME_getFirstElement( element ); child != NULL; child = xDataTOME_getNextElement( child ) ) {
        if( strcmp( child->name, "level" ) ) {
            smr_setReportError2( smr, smr_unknownID, 1, "invalid element '%s' in %s", child->name, element->name );
            goto err;
        }
        if( MCGIDI_target_heated_parseParticleLevel( smr, target, child, pop, mass_MeV, particleAliases ) ) goto err;
    }

    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_target_heated_parseParticleLevel( statusMessageReporting *smr, MCGIDI_target_heated *target, xDataTOM_element *element, MCGIDI_POP *parent,
    double mass_MeV, xDataTOM_element *particleAliases ) {

    int globalParticle = 0;
    char const *name, *level, *aliasValue;     /* Do not free any of these as they are owned by called routine. */
    double level_MeV = 0.;
    xDataTOM_element *alias;

    if( ( name = xDataTOM_getAttributesValueInElement( element, "name" ) ) == NULL ) {
        smr_setReportError2p( smr, smr_unknownID, 1, "particle missing name attribute" );
        goto err;
    }
    if( ( level = xDataTOM_getAttributesValueInElement( element, "energy" ) ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "particle '%s' level missing energy attribute", name );
        goto err;
    }
        /* Special case for 'c' labels. Correct mass is only needed for two-body. */
    if( level[0] != 'u' ) if( MCGIDI_misc_PQUStringToDoubleInUnitOf( smr, level, "MeV", &level_MeV ) ) goto err;
    for( alias = xDataTOME_getFirstElement( particleAliases ); alias != NULL; alias = xDataTOME_getNextElement( alias ) ) {
        if( ( aliasValue = xDataTOM_getAttributesValueInElement( alias, "value" ) ) == NULL ) {
            smr_setReportError2p( smr, smr_unknownID, 1, "particle missing name attribute" );
            goto err;
        }
        if( strcmp( aliasValue, name ) == 0 ) globalParticle = 1;
    }
    if( MCGIDI_POPs_addParticleIfNeeded( smr, &(target->pops), name, mass_MeV + level_MeV, level_MeV, parent, globalParticle ) == NULL ) goto err;

    return( MCGIDI_target_heated_parseParticleGammas( smr, target, element, name ) );

err:
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_target_heated_parseParticleGammas( statusMessageReporting *smr, MCGIDI_target_heated *target, xDataTOM_element *element, char const *name ) {

    int gammaCounts = 0;
    MCGIDI_POP *pop = MCGIDI_POPs_findParticle( &(target->pops), name );
    xDataTOM_element *child;
    MCGIDI_GammaBranching *gammas = NULL;
    char const *finalLevelString;
    double probability;

    for( child = xDataTOME_getFirstElement( element ); child != NULL; child = xDataTOME_getNextElement( child ) ) {
        if( strcmp( child->name, "gamma" ) ) {
            smr_setReportError2( smr, smr_unknownID, 1, "invalid element '%s' in %s", child->name, element->name );
            goto err;
        }
        gammaCounts++;
    }
    if( gammaCounts > 0 ) {
        if( ( gammas = (MCGIDI_GammaBranching *) smr_malloc2( smr, gammaCounts * sizeof( MCGIDI_GammaBranching), 0, "gammas" ) ) == NULL ) goto err;
        for( child = xDataTOME_getFirstElement( element ); child != NULL; child = xDataTOME_getNextElement( child ) ) {
            if( ( finalLevelString = xDataTOM_getAttributesValueInElement( child, "finalLevel" ) ) == NULL ) {
                smr_setReportError2p( smr, smr_unknownID, 1, "gamma missing 'finalLevel'" );
                goto err;
            }
            if( xDataTOME_convertAttributeToDouble( smr, child, "probability", &probability ) != 0 ) {
                smr_setReportError2p( smr, smr_unknownID, 1, "gamma missing 'probability' attribute" );
                goto err;
            }
        }
    }
    pop->numberOfGammaBranchs = gammaCounts;
    pop->gammas = gammas;

    return( 0 );

err:
    if( gammas != NULL ) smr_freeMemory( (void **) &gammas );
    return( 1 );
}
/*
************************************************************
*/
static int MCGIDI_target_heated_parseReaction( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_target_heated *target, 
        MCGIDI_POPs *pops, MCGIDI_reaction *reaction ) {

    if( MCGIDI_reaction_parseFromTOM( smr, element, target, pops, reaction ) ) goto err;
    target->numberOfReactions++;
    
    return( 0 );

err:
    smr_setReportError2( smr, smr_unknownID, 1, "%s\n", xDataTOM_getAttributesValueInElement( element, "outputChannel" ) );
    return( 1 );
}
/*
************************************************************
*/
int MCGIDI_target_heated_numberOfReactions( statusMessageReporting * /*smr*/, MCGIDI_target_heated *target ) {

    return( target->numberOfReactions );
}
/*
************************************************************
*/
int MCGIDI_target_heated_numberOfProductionReactions( statusMessageReporting * /*smr*/, MCGIDI_target_heated * /*target*/ ) {

    return( 0 );
}
/*
************************************************************
*/
MCGIDI_reaction *MCGIDI_target_heated_getReactionAtIndex( MCGIDI_target_heated *target, int index ) {

    if( ( index >= 0 ) && ( index < target->numberOfReactions ) ) return( &(target->reactions[index]) );
    return( NULL );
}
/*
************************************************************
*/
MCGIDI_reaction *MCGIDI_target_heated_getReactionAtIndex_smr( statusMessageReporting *smr, MCGIDI_target_heated *target, int index ) {

    MCGIDI_reaction *reaction = MCGIDI_target_heated_getReactionAtIndex( target, index );

    if( reaction == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "bad reaction index = %d for %s + %s", index, target->projectilePOP->name, target->targetPOP->name ); 
    }
    return( reaction );
}
#if 0
/*
************************************************************
*/
MCGIDI_channel *MCGIDI_target_heated_getProductionReactionAtIndex( MCGIDI_target_heated *target, int index ) {

    MCGIDI_channel *channel = NULL;

    if( ( index >= 0 ) && ( index < target->nProductionReactions ) ) channel = target->productionReactions[index];
    return( channel );
}
#endif
/*
************************************************************
*/
MCGIDI_POP *MCGIDI_target_heated_getPOPForProjectile( statusMessageReporting * /*smr*/, MCGIDI_target_heated *target ) {

    return( target->projectilePOP );
}
/*
************************************************************
*/
MCGIDI_POP *MCGIDI_target_heated_getPOPForTarget( statusMessageReporting * /*smr*/, MCGIDI_target_heated *target ) {

    return( target->targetPOP );
}
/*
************************************************************
*/
double MCGIDI_target_heated_getProjectileMass_MeV( statusMessageReporting * /*smr*/, MCGIDI_target_heated *target ) {

    return( MCGIDI_POP_getMass_MeV( target->projectilePOP ) );
}
/*
************************************************************
*/
double MCGIDI_target_heated_getTargetMass_MeV( statusMessageReporting * /*smr*/, MCGIDI_target_heated *target ) {

    return( MCGIDI_POP_getMass_MeV( target->targetPOP ) );
}
/*
************************************************************
*/
double MCGIDI_target_heated_getTotalCrossSectionAtE( statusMessageReporting *smr, MCGIDI_target_heated *target, 
        MCGIDI_quantitiesLookupModes &modes, bool sampling ) {

    double xsec;

    if( modes.getCrossSectionMode( ) == MCGIDI_quantityLookupMode_pointwise ) {
        double e_in = modes.getProjectileEnergy( );

        if( e_in < target->EMin ) e_in = target->EMin;
        if( e_in > target->EMax ) e_in = target->EMax;
        ptwXY_getValueAtX( target->crossSection, e_in, &xsec ); }
    else if( modes.getCrossSectionMode( ) == MCGIDI_quantityLookupMode_grouped ) {
        int index = modes.getGroupIndex( );
        double *xSecP;

        if( sampling ) {
            xSecP = ptwX_getPointAtIndex( target->crossSectionGroupedForSampling, index ); }
        else {
            xSecP = ptwX_getPointAtIndex( target->crossSectionGrouped, index );
        }

        if( xSecP != NULL ) {
            xsec = *xSecP; }
        else {
            xsec = 0.;
            smr_setReportError2( smr, smr_unknownID, 1, "Invalid cross section group index %d", index, (int) ptwX_length( target->crossSectionGrouped ) );
        } }
    else {
        xsec = 0.;
    }
    return( xsec );
}
/*
************************************************************
*/
double MCGIDI_target_heated_getIndexReactionCrossSectionAtE( statusMessageReporting *smr, MCGIDI_target_heated *target, int index, 
        MCGIDI_quantitiesLookupModes &modes, bool sampling ) {

    double xsec = 0.;
    MCGIDI_reaction *reaction = MCGIDI_target_heated_getReactionAtIndex_smr( smr, target, index );

    if( reaction != NULL ) xsec = MCGIDI_reaction_getCrossSectionAtE( smr, reaction, modes, sampling );
    return( xsec );
}
/*
************************************************************
*/
int MCGIDI_target_heated_sampleIndexReactionProductsAtE( statusMessageReporting *smr, MCGIDI_target_heated *target, int index, 
        MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo, MCGIDI_sampledProductsDatas *productDatas ) {

    MCGIDI_reaction *reaction = MCGIDI_target_heated_getReactionAtIndex_smr( smr, target, index );

    productDatas->numberOfProducts = 0;
    if( reaction == NULL ) return( -1 );
    return( MCGIDI_outputChannel_sampleProductsAtE( smr, &(reaction->outputChannel), modes, decaySamplingInfo, productDatas, NULL ) );
}
/*
************************************************************
*/
double MCGIDI_target_heated_getReactionsThreshold( statusMessageReporting * /*smr*/, MCGIDI_target_heated *target, int index ) {

    MCGIDI_reaction *reaction = MCGIDI_target_heated_getReactionAtIndex( target, index );

    if( reaction == NULL ) return( -1 );
    return( reaction->EMin );
}
/*
************************************************************
*/
int MCGIDI_target_heated_getReactionsDomain( statusMessageReporting * /*smr*/, MCGIDI_target_heated *target, int index, double *EMin, double *EMax ) {

    MCGIDI_reaction *reaction = MCGIDI_target_heated_getReactionAtIndex( target, index );

    if( reaction == NULL ) return( -1 );
    *EMin = reaction->EMin;
    *EMax = reaction->EMax;
    return( 0 );
}
/*
************************************************************
*/
double MCGIDI_target_heated_getIndexReactionFinalQ( statusMessageReporting *smr, MCGIDI_target_heated *target, int index, 
        MCGIDI_quantitiesLookupModes &modes ) {

    MCGIDI_reaction *reaction = MCGIDI_target_heated_getReactionAtIndex_smr( smr, target, index );

    if( reaction == NULL ) return( 0. );
    return( MCGIDI_reaction_getFinalQ( smr, reaction, modes ) );
}
/*
************************************************************
*/
std::map<int, enum MCGIDI_transportability> const *MCGIDI_target_heated_getUniqueProducts( statusMessageReporting * /*smr*/, MCGIDI_target_heated *target ) {

    return( target->transportabilities );
}
/*
************************************************************
*/
int MCGIDI_target_heated_recast( statusMessageReporting *smr, MCGIDI_target_heated *target, GIDI_settings &settings ) {

    int ir, projectilePoPID = target->projectilePOP->globalPoPsIndex;
    ptwXPoints *totalGroupedCrossSection = NULL;
    GIDI_settings_particle const *projectileSettings = settings.getParticle( projectilePoPID );
    nfu_status status_nf;

    if( projectileSettings == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "Settings missing for projectile %s", target->projectilePOP->name ); 
        return( 1 );
    }
    target->crossSectionGrouped = ptwX_free( target->crossSectionGrouped );
    target->crossSectionGroupedForSampling = ptwX_free( target->crossSectionGroupedForSampling );
    if( projectileSettings->isEnergyMode_grouped( ) ) {
        int64_t numberOfGroups = projectileSettings->getNumberOfGroups( );

        if( ( totalGroupedCrossSection = ptwX_createLine( numberOfGroups, numberOfGroups, 0, 0, &status_nf ) ) == NULL ) {
            smr_setReportError2( smr, smr_unknownID, 1, "totalGroupedCrossSection allocation failed: status_nf = %d, '%s'", 
                    status_nf, nfu_statusMessage( status_nf ) ); 
            goto err;
        }
    }

    for( ir = 0; ir < target->numberOfReactions; ++ir ) {
        if( MCGIDI_reaction_recast( smr, &(target->reactions[ir]), settings, projectileSettings, target->temperature_MeV, totalGroupedCrossSection ) != 0 ) goto err;
    }
    if( projectileSettings->isEnergyMode_grouped( ) ) {
        if( ( target->crossSectionGroupedForSampling = ptwX_clone( totalGroupedCrossSection, &status_nf ) ) == NULL ) {
            smr_setReportError2( smr, smr_unknownID, 1, "totalGroupedCrossSection allocation failed: status_nf = %d, '%s'", 
                    status_nf, nfu_statusMessage( status_nf ) ); 
            goto err;
        }
        for( ir = 0; ir < target->numberOfReactions; ++ir ) {
            int index = target->reactions[ir].thresholdGroupIndex;

            if( index > -1 ) {
                double xSec = target->reactions[ir].thresholdGroupedDeltaCrossSection +
                        ptwX_getPointAtIndex_Unsafely( target->crossSectionGroupedForSampling, index );

                ptwX_setPointAtIndex( target->crossSectionGroupedForSampling, index, xSec );
            }
        }
    }
    target->crossSectionGrouped = totalGroupedCrossSection;
    totalGroupedCrossSection = NULL;

    return( 0 );

err:
    ptwX_free( totalGroupedCrossSection );
    target->crossSectionGroupedForSampling = ptwX_free( target->crossSectionGroupedForSampling );
    return( 1 );
}

#if defined __cplusplus
}
#endif
