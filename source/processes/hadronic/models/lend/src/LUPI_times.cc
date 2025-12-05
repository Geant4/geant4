/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef _WIN32

#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include <iostream>

#include <LUPI.hpp>

namespace LUPI {

#define bufferSize 1024

/* *********************************************************************************************************//**
 * Default constructor for DeltaTime which sets all values to 0.0.
 ***********************************************************************************************************/

DeltaTime::DeltaTime( ) :
        m_CPU_time( 0.0 ),
        m_wallTime( 0.0 ),
        m_CPU_timeIncremental( 0.0 ),
        m_wallTimeIncremental( 0.0 ) {

}

/* *********************************************************************************************************//**
 * Constructor for DeltaTime which contains the CPU time *a_CPU_time* and wall time *a_wallTime*.
 ***********************************************************************************************************/

DeltaTime::DeltaTime( double a_CPU_time, double a_wallTime, double a_CPU_timeIncremental, double a_wallTimeIncremental ) :
        m_CPU_time( a_CPU_time ),
        m_wallTime( a_wallTime ),
        m_CPU_timeIncremental( a_CPU_timeIncremental ),
        m_wallTimeIncremental( a_wallTimeIncremental ) {

}

/* *********************************************************************************************************//**
 * Copy constructor for DeltaTime.
 ***********************************************************************************************************/

DeltaTime::DeltaTime( DeltaTime const &deltaTime ) :
        m_CPU_time( deltaTime.CPU_time( ) ),
        m_wallTime( deltaTime.wallTime( ) ),
        m_CPU_timeIncremental( deltaTime.CPU_timeIncremental( ) ),
        m_wallTimeIncremental( deltaTime.wallTimeIncremental( ) ) {

}

/* *********************************************************************************************************//**
 * Returns a string representation of *this*. The arugments *a_formatIncremental* and *a_format* specify, in sprintf style, the 
 * formats for the incremental and total times. If an arugment is an empty string (e.g., "") then its time is not included
 * in the return string. Each non-empty argument must contain two flags (e.g., "%.3f") for converting two doubles. If both formats 
 * are non-empty strings, then *a_sep* is inserted between them. An example of a format string is "total: CPU %8.3f, wall %8.3f".
 *
 * @param a_formatIncremental       [in]    Specifies the format in sprintf style for the incremental CPU and wall times.
 * @param a_formatTotal             [in]    Specifies the format in sprintf style for the total CPU and wall times.
 * @param a_sep                     [in]    Specifies the string that separates the total and incremental time strings.
 *
 * @return                                  A string representation of the delta times.
 ***********************************************************************************************************/

std::string DeltaTime::toString( std::string a_formatIncremental, std::string a_formatTotal, std::string a_sep ) {

    std::string deltaTimeStr;
    char buffer[bufferSize+1];

    if( a_formatIncremental != "" ) {
        snprintf( buffer, bufferSize, a_formatIncremental.c_str( ), m_CPU_timeIncremental, m_wallTimeIncremental );
        deltaTimeStr += buffer;
    }

    if( a_formatTotal != "" ) {
        if( deltaTimeStr != "" ) deltaTimeStr += a_sep;
        snprintf( buffer, bufferSize, a_formatTotal.c_str( ), m_CPU_time, m_wallTime );
        deltaTimeStr += buffer;
    }

    return( deltaTimeStr );
}

/* *********************************************************************************************************//**
 * Constructor.
 ***********************************************************************************************************/

Timer::Timer( ) {

    reset( );
}

/* *********************************************************************************************************//**
 * Returns a DeltaTime instance representing the time since *this* was created or *reset* was called.
 *  
 * @return                              A DeltaTime representing the time since *this* was created or *reset* was called.
 ***********************************************************************************************************/

DeltaTime Timer::deltaTime( ) {

    struct timeval wallTime;
    clock_t CPU_time = clock( );
    gettimeofday( &wallTime, 0 );

    double dWallTime = static_cast<double>( wallTime.tv_sec - m_wallTime.tv_sec ) 
            + 1e-6 * static_cast<double>( wallTime.tv_usec - m_wallTime.tv_usec );
    double dCPU_time = double( CPU_time - m_CPU_time ) / CLOCKS_PER_SEC;

    double dWallTimeIncremental = static_cast<double>( wallTime.tv_sec - m_wallTimeIncremental.tv_sec ) 
            + 1e-6 * static_cast<double>( wallTime.tv_usec - m_wallTimeIncremental.tv_usec );
    double dCPU_timeIncremental = double( CPU_time - m_CPU_timeIncremental ) / CLOCKS_PER_SEC;

    m_CPU_timeIncremental = CPU_time;
    m_wallTimeIncremental = wallTime;

    return( DeltaTime( dCPU_time, dWallTime, dCPU_timeIncremental, dWallTimeIncremental ) );
}

/* *********************************************************************************************************//**
 * Calls deltaTime and then reset. Returns the results of the call to deltaTime.
 *  
 * @return                              A DeltaTime representing the time since *this* was created or *reset* was called.
 ***********************************************************************************************************/

DeltaTime Timer::deltaTimeAndReset( ) {

    DeltaTime deltaTime1 = deltaTime( );
    reset( );

    return( deltaTime1 );
}

/* *********************************************************************************************************//**
 * Resets the internal times to the current time.
 ***********************************************************************************************************/

void Timer::reset( ) {

    m_CPU_time = clock( );
    gettimeofday( &m_wallTime, 0 );
    m_CPU_timeIncremental = m_CPU_time;
    m_wallTimeIncremental = m_wallTime;
}

}               // End of namespace LUPI.

#endif          // End of not _WIN32 defined.
