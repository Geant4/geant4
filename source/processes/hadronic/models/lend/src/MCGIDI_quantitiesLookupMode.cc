/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include "MCGIDI.h"

/*  ---- MCGIDI_quantitiesLookupModes ----  */
/*
=========================================================
*/
MCGIDI_quantitiesLookupModes::MCGIDI_quantitiesLookupModes( int projectilesPOPID ) {

    mProjectilesPOPID = projectilesPOPID;
    mProjectileEnergy = -1.;
    mGroupIndex = -1;
    mProjectileEnergyForGroupIndex = -1.;
    mTemperature = 0.;
    mCrossSectionMode = MCGIDI_quantityLookupMode_pointwise;
    mMultiplicityMode = MCGIDI_quantityLookupMode_pointwise;
}
/*
=========================================================
*/
MCGIDI_quantitiesLookupModes::~MCGIDI_quantitiesLookupModes( ) {

}
/*
=========================================================
*/
int MCGIDI_quantitiesLookupModes::setGroupIndex( GIDI_settings const &settings, bool encloseOutOfRange ) {

    GIDI_settings_particle const *particle = settings.getParticle( mProjectilesPOPID );
    if( particle == NULL ) throw 1;

    mGroupIndex = particle->getGroupIndexFromEnergy( mProjectileEnergy, encloseOutOfRange );
    if( mGroupIndex == -3 ) throw 1;

    mProjectileEnergyForGroupIndex = mProjectileEnergy;
    if( mGroupIndex < 0 ) mProjectileEnergyForGroupIndex = -1;
    return( mGroupIndex );
}
/*
=========================================================
*/
enum MCGIDI_quantityLookupMode MCGIDI_quantitiesLookupModes::getMode( std::string const &quantity ) const {

    if( quantity == std::string( "cross section" ) ) {
        return( mCrossSectionMode ); }
    else if( quantity == std::string( "multiplicity" ) ) {
        return( mMultiplicityMode ); }
    else {
        throw 1;
    }
}
/*
=========================================================
*/
std::vector<std::string> MCGIDI_quantitiesLookupModes::getListOfLookupQuanities( ) const {

    std::vector<std::string> quanities;

    quanities.push_back( std::string( "cross section" ) );
    quanities.push_back( std::string( "multiplicity" ) );

    return( quanities );
}
/*
=========================================================
*/
void MCGIDI_quantitiesLookupModes::setMode( std::string const &quantity, enum MCGIDI_quantityLookupMode mode ) {

    if( quantity == std::string( "cross section" ) ) {
        mCrossSectionMode = mode; }
    else if( quantity == std::string( "multiplicity" ) ) {
        mMultiplicityMode = mode; }
    else {
        throw 1;
    }
}
/*
=========================================================
*/
void MCGIDI_quantitiesLookupModes::setModeAll( enum MCGIDI_quantityLookupMode mode ) {

    mCrossSectionMode = mode;
    mMultiplicityMode = mode;
}
