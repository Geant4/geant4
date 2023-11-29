/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include "MCGIDI.h"

static char versionStr[64] = "";

/*
========================================================================
*/
const char *MCGIDI_version( void ) {

    if( versionStr[0] == 0 ) snprintf( versionStr, sizeof versionStr, "MCGIDI version %d.%d.%d", MCGIDI_VERSION_MAJOR, MCGIDI_VERSION_MINOR, MCGIDI_VERSION_PATCHLEVEL );
    return( versionStr );
}
/*
========================================================================
*/
int MCGIDI_versionMajor( void ) { return( MCGIDI_VERSION_MAJOR ); }
int MCGIDI_versionMinor( void ) { return( MCGIDI_VERSION_MINOR ); }
int MCGIDI_versionPatchLevel( void ) { return( MCGIDI_VERSION_PATCHLEVEL ); }
