/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>

#include "ptwX.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
void ptwX_simpleWrite( ptwXPoints const *ptwX, FILE *f, char const *format ) {

    int64_t i1;
    double *p1 = ptwX->points;

    for( i1 = 0; i1 < ptwX->length; ++i1, ++p1 ) fprintf( f, format, *p1 );
}
/*
************************************************************
*/
void ptwX_simplePrint( ptwXPoints const *ptwX, char const *format ) {

    ptwX_simpleWrite( ptwX, stdout, format );
}

#if defined __cplusplus
}
#endif
