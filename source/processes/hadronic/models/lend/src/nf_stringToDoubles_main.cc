/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nf_utilities.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif
/*
char str[] = "1e-5 20.43634 2e-5 20.43634 5e-5 20.43634 1e-4 20.43633 2e-4 20.43633 5e-4 20.43633 1e-3 20.43633 2e-3 20.43633 5e-3 20.43633 " \
"1e-2 20.43633 0.0253 20.43633 5e-2 20.43633 0.1 20.43632 0.2 20.43631 0.5 20.43627 1 20.4362 2 20.43606 5 20.43566 10 20.43499 20" \
" 20.43364 50 20.4296 1e2 20.42288 2e2 20.40944 5e2 20.36926 1e3 20.30269 2e3 20.17105 4e3 19.91352 6e3 19.66341 8e3 19.4204 1e4 " \
"19.18418 1.5e4 18.62156 2e4 18.09567 4e4 16.29554 6e4 14.86728 8e4 13.70583 1e5 12.74229 1.5e5 10.9234 2e5 9.643178 3e5 7.951949 " \
"4e5 6.876412 5e5 6.125445 6e5 5.566879 7e5 5.13201 8e5 4.78157 9e5 4.491471 1e6 4.246104 1.2e6 3.850454 1.4e6 3.541748 1.6e6 3.291314 "
" 1.8e6 3.082187 2e6 2.903645 2.2e6 2.748543 2.4e6 2.611918 2.6e6 2.490197 2.8e6 2.380736 3e6 2.281521 3.2e6 2.190993 3.4e6 2.107917 3.6e6 " \
"2.031301 3.8e6 1.960334 4e6 1.894349 4.2e6 1.832787 4.4e6 1.775177 4.6e6 1.721118 4.8e6 1.670264 5e6 1.622318 5.5e6 1.513553 6e6 1.418157 " \
"6.5e6 1.333709 7e6 1.258367 7.5e6 1.190697 8e6 1.129564 8.5e6 1.074052 9e6 1.023415 9.5e6 0.9770347 1e7 0.9343974 1.05e7 0.8950685 1.1e7 " \
"0.8586796 1.15e7 0.8249154 1.2e7 0.7935044 1.25e7 0.7642113 1.3e7 0.7368313 1.35e7 0.7111848 1.4e7 0.6871141 1.45e7 0.6644799 1.5e7 0.6431586 " \
"1.55e7 0.6230401 1.6e7 0.6040262 1.65e7 0.586029 1.7e7 0.5689692 1.75e7 0.5527757 1.8e7 0.537384 1.85e7 0.5227359 1.9e7 0.5087783 1.95e7 0.495463 " \
"2e7 0.4827462 ";
*/

/*
========================================================================
*/
/*
int main( int argc, char **argv ) {

    int64_t i1, numberConverted;
    double *doublePtr;
    nfu_status status;
    char *endCharacter;

    status = nfu_stringToListOfDoubles( str, &numberConverted, &doublePtr, &endCharacter );
    if( doublePtr != NULL ) {
        for( i1 = 0; i1 < numberConverted; i1++ ) printf( "%6d %14.7e\n", (int) i1, doublePtr[i1] );
        nfu_free( doublePtr );
    }
    printf( "No converted = <%s>\n", endCharacter );
    printf( "%8d  %d\n", (int) numberConverted, status );
    printf( "%s\n", nfu_statusMessage( status ) );
    exit( EXIT_SUCCESS );
}
*/

#if defined __cplusplus
}
#endif
