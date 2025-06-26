/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>

#include "ptwX.h"

/*
************************************************************
*/
nfu_status ptwX_simpleWrite( statusMessageReporting *smr, ptwXPoints const *ptwX, FILE *f, char const *format ) {

    int64_t i1;
    double *p1 = ptwX->points;

    for( i1 = 0; i1 < ptwX->length; ++i1, ++p1 ) fprintf( f, format, *p1 );
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_simplePrint( statusMessageReporting *smr, ptwXPoints const *ptwX, char const *format ) {

    return( ptwX_simpleWrite( smr, ptwX, stdout, format ) );
}
