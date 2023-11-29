/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <iostream>
#include <stdlib.h>

#include "GIDI_settings.hh"

/*
=========================================================
*/
GIDI_settings_group::GIDI_settings_group( std::string const &label, int size1 ) {

    initialize( label, size1, size1, NULL );
}
/*
=========================================================
*/
GIDI_settings_group::GIDI_settings_group( std::string const &label, int length, double const *boundaries ) {

    initialize( label, length, length, boundaries );
}
/*
=========================================================
*/
GIDI_settings_group::GIDI_settings_group( std::string const &label, std::vector<double> const &boundaries) {

    int size1 = (int) boundaries.size( );

    initialize( label, size1, size1, &(boundaries[0]) );
}
/*
=========================================================
*/
GIDI_settings_group::GIDI_settings_group( GIDI_settings_group const &group ) {

    initialize( group.mLabel, group.size( ), group.size( ), &(group.mBoundaries[0]) );
}
/*
=========================================================
*/
void GIDI_settings_group::initialize( std::string const &label, int size1, int length, double const *boundaries ) {

    int i1;

    mLabel = label;
    if( size1 < length ) size1 = length;
    if( size1 < 0 ) size1 = 0;
    mBoundaries.resize( size1, 0 );
    for( i1 = 0; i1 < length; ++i1 ) mBoundaries[i1] = boundaries[i1];
}
/*
=========================================================
*/
GIDI_settings_group& GIDI_settings_group::operator=( const GIDI_settings_group &group ) {
  if ( this != &group ) {
    initialize( group.mLabel, group.size(), group.size(), &(group.mBoundaries[0]) );
  }
  return *this;
}
/*
=========================================================
*/
GIDI_settings_group::~GIDI_settings_group( ) {

}
/*
=========================================================
*/
int GIDI_settings_group::getGroupIndexFromEnergy( double energy, bool encloseOutOfRange ) const {

    int iMin = 0, iMid, iMax = (int) mBoundaries.size( ), iMaxM1 = iMax - 1;

    if( iMax == 0 ) return( -3 );
    if( energy < mBoundaries[0] ) {
        if( encloseOutOfRange ) return( 0 );
        return( -2 );
    }
    if( energy > mBoundaries[iMaxM1] ) {
        if( encloseOutOfRange ) return( iMax - 2 );
        return( -1 );
    }
    while( 1 ) { // Loop checking, 11.06.2015, T. Koi
        iMid = ( iMin + iMax ) >> 1;
        if( iMid == iMin ) break;
        if( energy < mBoundaries[iMid] ) {
            iMax = iMid; }
        else {
            iMin = iMid;
        }
    }
    if( iMin == iMaxM1 ) iMin--;
    return( iMin );
}
/*
=========================================================
*/
void GIDI_settings_group::print( bool outline, int valuesPerLine ) const {

    int nbs = size( );
    char buffer[128];

    std::cout << "GROUP: label = '" << mLabel << "': length = " << nbs << std::endl;
    if( outline ) return;
    for( int ib = 0; ib < nbs; ib++ ) {
        snprintf( buffer, sizeof buffer, "%16.8e", mBoundaries[ib] );
        std::cout << buffer;
        if( ( ( ib + 1 ) % valuesPerLine ) == 0 ) std::cout << std::endl;
    }
    if( nbs % valuesPerLine ) std::cout << std::endl;
}

#if 0
/*  ---- GIDI_settings_groups_from_bdfls ----  */
/*
=========================================================
*/
GIDI_settings_groups_from_bdfls::GIDI_settings_groups_from_bdfls( std::string const &fileName ) {

    initialize( fileName.c_str( ) );
}
/*
=========================================================
*/
GIDI_settings_groups_from_bdfls::GIDI_settings_groups_from_bdfls( char const *fileName ) {

    initialize( fileName );
}
/*
=========================================================
*/
GIDI_settings_groups_from_bdfls::GIDI_settings_groups_from_bdfls( cbdfls_file const *bdfls ) {

    initialize2( bdfls );
}
/*
=========================================================
*/
void GIDI_settings_groups_from_bdfls::initialize( char const *fileName ) {

    cbdfls_file *bdfls;
    cbdflsErrors Error;

    if( ( bdfls = cbdflsOpen( fileName, &Error ) ) == NULL ) throw Error;
    initialize2( bdfls );
    cbdflsRelease( bdfls );
}
/*
=========================================================
*/
void GIDI_settings_groups_from_bdfls::initialize2( cbdfls_file const *bdfls ) {

    int ng, ngbs, *gids;
    double *boundaries;
    std::string label( "" );
    char cLabel[100];

    ng = cbdflsGIDs( (cbdfls_file *) bdfls, &gids );
    for( int ig = 0; ig < ng; ++ig ) {
        ngbs = cbdflsGetGroup( (cbdfls_file *) bdfls, gids[ig], &boundaries );
        snprintf( cLabel, sizeof xLabel, "LLNL_gid_%.3d", gids[ig] );
        label = cLabel;
        mGroups.push_back( GIDI_settings_group( label, ngbs, boundaries ) );
    }
}
/*
=========================================================
*/
GIDI_settings_groups_from_bdfls::~GIDI_settings_groups_from_bdfls( ) {

}
/*
=========================================================
*/
GIDI_settings_group GIDI_settings_groups_from_bdfls::getViaGID( int gid ) const {

    std::string label( "" );
    char cLabel[100];

    snprintf( cLabel, sizeof cLabel, "LLNL_gid_%.3d", gid );
    label = cLabel;
    for( int ig = 0; ig < (int) mGroups.size( ); ++ig ) {
        if( mGroups[ig].isLabel( label ) ) return( mGroups[ig] );
    }
    throw 1;
}
/*
=========================================================
*/
std::vector<std::string> GIDI_settings_groups_from_bdfls::getLabels( void ) const {

    int size = (int) mGroups.size( );
    std::vector<std::string> labels( size );

    for( int if1 = 0; if1 < size; ++if1 ) labels[if1] = mGroups[if1].getLabel( );
    return( labels );
}
/*
=========================================================
*/
std::vector<int> GIDI_settings_groups_from_bdfls::getGIDs( void ) const {

    int size = (int) mGroups.size( );
    std::vector<int> fids( size );
    char *e;

    for( int if1 = 0; if1 < size; ++if1 ) {
        fids[if1] = (int) strtol( &(mGroups[if1].getLabel( ).c_str( )[9]), &e, 10 );
    }
    return( fids );
}
/*
=========================================================
*/
void GIDI_settings_groups_from_bdfls::print( bool outline, int valuesPerLine ) const {

    int ngs = (int) mGroups.size( );

    std::cout << "BDFLS GROUPs: number of groups = " << ngs << std::endl;
    for( int if1 = 0; if1 < ngs ; ++if1 ) mGroups[if1].print( outline, valuesPerLine );
}
#endif
