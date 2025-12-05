/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <iostream>
#include <iomanip>

#include <LUPI.hpp>

namespace LUPI {

/* *********************************************************************************************************//**
 * Constructor for a StatusMessageReporting instance.
 ***********************************************************************************************************/

StatusMessageReporting::StatusMessageReporting( ) {

    int status = smr_initialize( &m_smr, smr_status_Ok );

    if( status != 0 ) throw( "StatusMessageReporting::StatusMessageReporting: Oops." ); // Currently, this should never happend.
}

/* *********************************************************************************************************//**
 * Destructor for a StatusMessageReporting instance.
 ***********************************************************************************************************/

StatusMessageReporting::~StatusMessageReporting( ) {

    smr_release( &m_smr );
}

/* *********************************************************************************************************//**
 * Returns the first *a_reports* reports from *m_smr* with *a_prefix* appended to the beginning of the returned string.
 *
 * @param a_prefix          [in]    A string added to the beginning of the message.
 * @param a_report          [in]    The maximum number of reports to include in the message.
 * @param a_clear           [in]    If *true*, calls the **clear()** method after the message is constructed.
 ***********************************************************************************************************/

std::string StatusMessageReporting::constructMessage( std::string a_prefix, int a_reports, bool a_clear ) {

    std::string sep( "" );
    std::string message( a_prefix );
    statusMessageReport const *report;

    if( a_prefix == "" ) sep = "\n";

    for( report = smr_firstReport( &m_smr ); report != NULL; report = smr_nextReport( report ), --a_reports ) {
        if( a_reports == 0 ) break;

        char *reportMessage = smr_copyMessage( report );
        if( reportMessage != nullptr ) {
            message += sep;
            message += reportMessage;
            free( reportMessage );
            sep = "\n";
        }
    }
    if( a_clear ) clear( );

    return( message );
}

/* *********************************************************************************************************//**
 * Returns the first *a_reports* reports from *m_smr* with *a_prefix* appended to the beginning of the returned string.
 *
 * @param a_prefix          [in]    A string added to the beginning of the message.
 * @param a_report          [in]    The maximum number of reports to include in the message.
 * @param a_clear           [in]    If *true*, calls the **clear()** method after the message is constructed.
 ***********************************************************************************************************/

std::string StatusMessageReporting::constructFullMessage( std::string const &a_prefix, int a_reports, bool a_clear ) {

    std::string message( a_prefix );
    statusMessageReport const *report;

    for( report = smr_firstReport( &m_smr ); report != NULL; report = smr_nextReport( report ), --a_reports ) {
        if( a_reports == 0 ) break;

        char *reportMessage = smr_copyFullMessage( report );
        if( reportMessage != nullptr ) {
            message += '\n';
            message += reportMessage;
            free( reportMessage );
        }
    }
    if( a_clear ) clear( );

    return( message );
}

}               // End of namespace LUPI.
