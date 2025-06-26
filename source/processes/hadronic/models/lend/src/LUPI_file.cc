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

#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#include <filesystem>
char *realpath( char const *a_path, char *a_resolved ) {

    char resolvedPath[LUPI_PATH_MAX+1], *p1 = nullptr;

    DWORD length = GetFullPathName( a_path, LUPI_PATH_MAX, resolvedPath, nullptr );
    // MSVC requires explicitly casting malloc result
    if( ( p1 = (char *)malloc( length + 1 ) ) == nullptr ) return( nullptr );
    strcpy( p1, resolvedPath );
    if( length == 0 ) return( nullptr );
    return( p1 );
}

std::string dirname( char const *a_path ) {
    std::filesystem::path filePath( a_path );
    return filePath.parent_path().string();
}

std::string basename( char const *a_path ) {
    // don't strip the extension for compatibility with <unistd> version
    std::filesystem::path filePath( a_path );
    return filePath.filename().string();
}
#else
// FIXME: once all users are on C++17 or later, switch to using std::filesystem for all systems
#include <unistd.h>
#include <libgen.h>
#endif

namespace LUPI {

namespace FileInfo {

/* *********************************************************************************************************//**
 * This function takes a file path and returns its real path. On a Unix system, the system function realpath is called.
 *      
 * @param a_path        [in]    The path whose real path is to be determined.
 *      
 * @return                      The real path.
 ***********************************************************************************************************/
        
std::string realPath( std::string const &a_path ) {    
        
    char *p1 = realpath( a_path.c_str( ), nullptr );
        
    if( p1 == nullptr ) {
        std::string errMsg( "realPath: file does not exist: " );
        throw Exception( errMsg + a_path );
    } 
    std::string basePath( p1 );
    free( p1 );
    return( basePath );
}

/* *********************************************************************************************************//**
 * Returns the base name of a path.
 *
 * @param a_path            [in]    The path whose base name is returned.
 ***********************************************************************************************************/

std::string _basename( std::string const &a_path ) {

    char *path = new char[a_path.size( ) + 1];
    strcpy( path, a_path.c_str( ) );
    std::string basename1( basename( path ) );

    delete[] path;

    return( basename1 );
}

/* *********************************************************************************************************//**
 * Same as **basename** but removes, if a *period* (i.e., ".") exists in the string, the last "." and all following characters.
 *
 * @param a_path            [in]    The path whose base name is returned without its extension.
 ***********************************************************************************************************/

std::string basenameWithoutExtension( std::string const &a_path ) {

    std::size_t found = a_path.rfind( '.' );

    return( a_path.substr( 0, found ) );
}

/* *********************************************************************************************************//**
 * Returns the directory name of a path. This is, it removes the base name.
 *
 * @param a_path            [in]    The path whose base name is returned.
 ***********************************************************************************************************/

std::string _dirname( std::string const &a_path ) {

    char *path = new char[a_path.size( ) + 1];
    strcpy( path, a_path.c_str( ) );
    std::string dirname1( dirname( (char *) path ) );

    delete[] path;

    return( dirname1 );
}

/* *********************************************************************************************************//**
 * Returns *true* if the path exists and *false* otherwise.
 *
 * @param a_path            [in]    The path that is checked for existence.
 ***********************************************************************************************************/

bool exists( std::string const &a_path ) {

#ifdef _WIN32
        return std::filesystem::exists( std::filesystem::path( a_path ) );
#else
        return( access( a_path.c_str( ), F_OK ) == 0 );
#endif
}

/* *********************************************************************************************************//**
 * Returns *true* if path is a direction that exists and *false* otherwise.
 *
 * @param a_path            [in]    The path that is checked for existence and is it a directory.
 *
 * @return                          Returns *true* if the path exists (e.g., created if it does not exists) and *false* otherwise.
 ***********************************************************************************************************/

bool isDirectory( std::string const &a_path ) {

    try {
        FileStat fileStat( a_path );
        return( fileStat.isDirectory( ) ); }
    catch (...) {
    }

    return( false );
}

/* *********************************************************************************************************//**
 * Adds all needed directories to complete *a_path*.
 *
 * @param a_path            [in]    The path that is checked for existence.
 *
 * @return                          Returns *true* if the path exists (e.g., created if it does not exists) and *false* otherwise.
 ***********************************************************************************************************/

bool createDirectories( std::string const &a_path ) {

    if( isDirectory( a_path ) ) return( true );
    if( ( a_path == LUPI_FILE_SEPARATOR ) || ( a_path == "." ) || ( a_path == "" ) ) return( true );

    std::string dirname1( _dirname( a_path ) );
    if( createDirectories( dirname1 ) ) {
#ifdef _WIN32
        int status = _mkdir( a_path.c_str( ) );
#else
        int status = mkdir( a_path.c_str( ), S_IRWXU | S_IRWXG | S_IRWXG );
#endif
        if( status == 0 ) return( true );
        switch( errno ) {
        case EEXIST :
            return( true );
        default :
            return( false );
        }
    }

    return( false );
}

/* *********************************************************************************************************//**
 * Calls the C stat function and stores its information.
 *
 * @param a_path            [in]    The path (e.g., file, directory) whose stat is determined.
 ***********************************************************************************************************/

FileStat::FileStat( std::string const &a_path ) :
        m_path( a_path ) {

    int error = stat( a_path.c_str( ), &m_stat );

    if( error != 0 ) {
        switch( error ) {
        case EACCES :
            throw Exception( "FileStat::FileStat: Permission denied for file '" + a_path + "'.." );
        case EIO :
            throw Exception( "FileStat::FileStat: An error occurred while stat-ing file '" + a_path + "'.." );
        case ELOOP :
            throw Exception( "FileStat::FileStat: A loop exists in symbolic links for file '" + a_path + "'.." );
        case ENAMETOOLONG :
            throw Exception( "FileStat::FileStat: Path name too long '" + a_path + "'." );
        case ENOENT :
            throw Exception( "FileStat::FileStat: No such path '" + a_path + "'." );
        case ENOTDIR :
            throw Exception( "FileStat::FileStat: A component of the path prefix is not a directory '" + a_path + "'." );
        case EOVERFLOW :
            throw Exception( "FileStat::FileStat: File too big: '" + a_path + "'." );
        default :
            throw Exception( "FileStat::FileStat: Unknown error from C function 'stat' for file '" + a_path + "'." );
        }
    }
}

}               // End of namespace FileInfo.

}               // End of namespace LUPI.
