/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <iostream>
#include <stdlib.h>

#include "GIDI_settings.hh"

/*  ---- GIDI_settings_flux_order ----  */
/*
=========================================================
*/
GIDI_settings_flux_order::GIDI_settings_flux_order( int order ) {

    if( order < 0 ) throw 1;
    mOrder = order;
}
/*
=========================================================
*/
GIDI_settings_flux_order::GIDI_settings_flux_order( int order, int length, double const *energies, double const *fluxes ) {

    initialize( order, length, energies, fluxes );
}
/*
=========================================================
*/
GIDI_settings_flux_order::GIDI_settings_flux_order( int order, std::vector<double> const &energies, std::vector<double> const &fluxes ) {

    int length = (int) energies.size( );

    if( length != (int) fluxes.size( ) ) throw 1;
    initialize( order, length, &(energies[0]), &(fluxes[0]) );
}
/*
=========================================================
*/
GIDI_settings_flux_order::GIDI_settings_flux_order( GIDI_settings_flux_order const &fluxOrder ) {

    initialize( fluxOrder.mOrder, fluxOrder.size( ), &(fluxOrder.mEnergies[0]), &(fluxOrder.mFluxes[0]) );
}
/*
=========================================================
*/
void GIDI_settings_flux_order::initialize( int order, int length, double const *energies, double const *fluxes ) {

    if( order < 0 ) throw 1;
    mOrder = order;
    mEnergies.resize( length, 0 );
    mFluxes.resize( length, 0 );
    for( int i1 = 0; i1 < length; ++i1 ) mEnergies[i1] = energies[i1];
    for( int i1 = 0; i1 < length; ++i1 ) mFluxes[i1] = fluxes[i1];
}
/*
=========================================================
*/
GIDI_settings_flux_order& GIDI_settings_flux_order::operator=( const GIDI_settings_flux_order &fluxOrder ) {
  if ( this != &fluxOrder ) {
    initialize( fluxOrder.mOrder, fluxOrder.size(), &(fluxOrder.mEnergies[0]), &(fluxOrder.mFluxes[0]) );
  }
  return *this;
}
/*
=========================================================
*/
GIDI_settings_flux_order::~GIDI_settings_flux_order( ) {

}
/*
=========================================================
*/
void GIDI_settings_flux_order::print( int valuesPerLine ) const {

    int nE = (int) mEnergies.size( );
    bool printIndent = true;
    char buffer[2 * 128];

    std::cout << "    ORDER: " << mOrder << std::endl;
    for( int iE = 0; iE < nE; ++iE ) {
        if( printIndent ) std::cout << "    ";
        printIndent = false;
        snprintf( buffer, sizeof buffer, "   %15.8e %15.8e", mEnergies[iE], mFluxes[iE] );
        std::cout << buffer;
        if( ( ( iE + 1 ) % valuesPerLine ) == 0 ) {
            std::cout << std::endl;
            printIndent = true;
        }
    }
    if( nE % valuesPerLine ) std::cout << std::endl;
}

/*  ---- GIDI_settings_flux ----  */
/*
=========================================================
*/
GIDI_settings_flux::GIDI_settings_flux( std::string const &label, double temperature ) {

    mLabel = label;
    mTemperature = temperature;
}
/*
=========================================================
*/
GIDI_settings_flux::GIDI_settings_flux( char const *label, double temperature ) {

    mLabel = label;
    mTemperature = temperature;
}
/*
=========================================================
*/
GIDI_settings_flux::GIDI_settings_flux( GIDI_settings_flux const &flux ) {

    mLabel = flux.getLabel( );
    mTemperature = flux.mTemperature;
    for( std::vector<GIDI_settings_flux_order>::const_iterator iter = flux.mFluxOrders.begin( ); iter < flux.mFluxOrders.end( ); ++iter ) addFluxOrder( *iter );
}
/*
=========================================================
*/
GIDI_settings_flux& GIDI_settings_flux::operator=( const GIDI_settings_flux &flux ) {
  if ( this != &flux ) {
    mLabel = flux.getLabel();
    mTemperature = flux.mTemperature;
    for( std::vector<GIDI_settings_flux_order>::const_iterator iter = flux.mFluxOrders.begin( ); iter < flux.mFluxOrders.end( ); ++iter ) addFluxOrder( *iter );
  }
  return *this;
}
/*
=========================================================
*/
GIDI_settings_flux::~GIDI_settings_flux( ) {

}
/*
=========================================================
*/
GIDI_settings_flux_order const *GIDI_settings_flux::operator[]( int index ) const {

    return( &(mFluxOrders[index]) );
}
/*
=========================================================
*/
void GIDI_settings_flux::addFluxOrder( GIDI_settings_flux_order const &fluxOrder ) {
/*
*   Orders can only be added in sequence (e.g., 0 first, then 1, ...).
*/
    int order = fluxOrder.getOrder( );

    if( order > (int) mFluxOrders.size( ) ) throw 1;
    mFluxOrders.push_back( fluxOrder );
}
/*
=========================================================
*/
void GIDI_settings_flux::print( bool outline, int valuesPerLine ) const {

    std::cout << "FLUX: label = '" << mLabel << "': maximum order = " << ( size( ) + 1 ) << std::endl;
    if( outline ) return;
    for( std::vector<GIDI_settings_flux_order>::const_iterator iter = mFluxOrders.begin( ); iter < mFluxOrders.end( ); ++iter ) iter->print( valuesPerLine );
}

#if 0
/*  ---- GIDI_settings_fluxes_from_bdfls ----  */
/*
=========================================================
*/
GIDI_settings_fluxes_from_bdfls::GIDI_settings_fluxes_from_bdfls( std::string const &fileName, double temperature_MeV = 0 ) {

    initialize( fileName.c_str( ), temperature_MeV );
}
/*
=========================================================
*/
GIDI_settings_fluxes_from_bdfls::GIDI_settings_fluxes_from_bdfls( char const *fileName, double temperature_MeV = 0 ) {

    initialize( fileName, temperature_MeV );
}
/*
=========================================================
*/
GIDI_settings_fluxes_from_bdfls::GIDI_settings_fluxes_from_bdfls( cbdfls_file const *bdfls, double temperature_MeV = 0 ) {

    initialize2( bdfls, temperature_MeV );
}
/*
=========================================================
*/
void GIDI_settings_fluxes_from_bdfls::initialize( char const *fileName, double temperature_MeV ) {

    cbdfls_file *bdfls;
    cbdflsErrors Error;

    if( ( bdfls = cbdflsOpen( fileName, &Error ) ) == NULL ) throw Error;
    initialize2( bdfls, temperature_MeV );
    cbdflsRelease( bdfls );
}
/*
=========================================================
*/
void GIDI_settings_fluxes_from_bdfls::initialize2( cbdfls_file const *bdfls, double temperature_MeV ) {

    int nf, length, *fids, order;
    double *energies, *fluxes;
    char label[100];

    nf = cbdflsFIDs( (cbdfls_file *) bdfls, &fids );
    for( int if1 = 0; if1 < nf; ++if1 ) {
        snprintf( label, sizeof label, "LLNL_fid_%.3d", fids[if1] );
        GIDI_settings_flux flux = GIDI_settings_flux( label, temperature_MeV );
        order = cbdflsGetFluxOrder( (cbdfls_file *) bdfls, fids[if1] );
        for( int io = 0; io <= order; ++io ) {
            length = cbdflsGetFlux( (cbdfls_file *) bdfls, fids[if1], io, &energies, &fluxes );
            GIDI_settings_flux_order flux_order = GIDI_settings_flux_order( io, length, energies, fluxes );
            flux.addFluxOrder( flux_order );
        }
        mFluxes.push_back( flux );
    }
    return;
}
/*
=========================================================
*/
GIDI_settings_fluxes_from_bdfls::~GIDI_settings_fluxes_from_bdfls( ) {

}
/*
=========================================================
*/
GIDI_settings_flux GIDI_settings_fluxes_from_bdfls::getViaFID( int fid ) {

    char label[100];

    snprintf( label, sizeof label, "LLNL_fid_%.3d", fid );
    for( int if1 = 0; if1 < (int) mFluxes.size( ); ++if1 ) {
        if( mFluxes[if1].isLabel( label ) ) return( mFluxes[if1] );
    }
    throw 1;
}
/*
=========================================================
*/
std::vector<std::string> GIDI_settings_fluxes_from_bdfls::getLabels( void ) {

    int size = (int) mFluxes.size( );
    std::vector<std::string> labels( size );

    for( int if1 = 0; if1 < size; ++if1 ) labels[if1] = mFluxes[if1].getLabel( );
    return( labels );
}
/*
=========================================================
*/
std::vector<int> GIDI_settings_fluxes_from_bdfls::getFIDs( void ) {

    int size = (int) mFluxes.size( );
    std::vector<int> fids( size );
    char *e;

    for( int if1 = 0; if1 < size; ++if1 ) {
        fids[if1] = (int) strtol( &(mFluxes[if1].getLabel( ).c_str( )[9]), &e, 10 );
    }
    return( fids );
}
/*
=========================================================
*/
void GIDI_settings_fluxes_from_bdfls::print( bool outline, int valuesPerLine ) {

    int nfs = (int) mFluxes.size( );

    std::cout << "BDFLS FLUXes: number of fluxes = " << nfs << std::endl;
    for( int if1 = 0; if1 < nfs ; ++if1 ) mFluxes[if1].print( outline, valuesPerLine );
}
#endif
