/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>

#include "xDataTOM_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
int xDataTOM_LegendreSeries_initialize( statusMessageReporting *smr, xDataTOM_LegendreSeries *LegendreSeries, int index, int length, double value ) {

    LegendreSeries->LegendreSeries = NULL;
    LegendreSeries->index = index;
    LegendreSeries->length = length;
    LegendreSeries->value = value;
    if( ( LegendreSeries->LegendreSeries = (double *) smr_malloc2( smr, length * sizeof( double ), 0, "LegendreSeries->LegendreSeries" ) ) == NULL ) return( 1 );

    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_LegendreSeries_release( xDataTOM_LegendreSeries *LegendreSeries ) {

    if( LegendreSeries == NULL ) return( 0 );
    smr_freeMemory( (void **) &(LegendreSeries->LegendreSeries) );
    return( 0 );
}

#if defined __cplusplus
}
#endif
