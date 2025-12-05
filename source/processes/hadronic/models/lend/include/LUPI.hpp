/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef LUPI_hpp_included
#define LUPI_hpp_included 1

#include <sys/stat.h>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <stdexcept>
#include <iostream>
#ifndef _WIN32
    #include <time.h>
    #include <sys/time.h>
#endif

#include <LUPI_defines.hpp>
#include <statusMessageReporting.h>

#define LUPI_XML_verionEncoding "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"

#ifndef LUPI_PATH_MAX
#define LUPI_PATH_MAX ( 4 * 4096 )
#endif

#if defined (GIDIP_HAVE_COMPILER_FLOATING_POINT_EXCEPTIONS)
void LUPI_FPE_enable( char const *a_file, int a_line );
void LUPI_FPE_disable_and_clear( char const *a_file, int a_line );
void LUPI_FPE_test( char const *a_file, int a_line );
#endif

namespace LUPI {

#ifdef _WIN32
#define LUPI_FILE_SEPARATOR   "\\"
#else
#define LUPI_FILE_SEPARATOR   "/"
#endif

#define GNDS_formatVersion_1_10Chars "1.10"
#define GNDS_formatVersion_2_0Chars "2.0"
#define GNDS_formatVersion_2_0_LLNL_4Chars "2.0.LLNL_4"
#define GNDS_formatVersion_2_1Chars "2.1"

void deprecatedFunction( std::string const &a_functionName, std::string const &a_replacementName, std::string const &a_asOf );

/*
============================================================
====================== FormatVersion =======================
============================================================
*/
class FormatVersion {

    private:
        std::string m_format;               /**< The GNDS format version. */
        int m_major;                        /**< The GNDS format major value as an integer. */
        int m_minor;                        /**< The GNDS format minor value as an integer. */
        std::string m_patch;                /**< The GNDS format patch string. This will be an empty string except for unofficial formats. */

    public:
        FormatVersion( );
        FormatVersion( std::string const &a_formatVersion );
        FormatVersion( FormatVersion const &a_formatVersion );
        FormatVersion &operator=( FormatVersion const &a_rhs );

        std::string const &format( ) const { return( m_format ); }
        int major( ) const { return( m_major ); }
        int minor( ) const { return( m_minor ); }
        std::string const &patch( ) const { return( m_patch ); }

        bool setFormat( std::string const &a_formatVersion );
        bool supported( ) const ;
};

/*
============================================================
========================= Exception ========================
============================================================
*/
class Exception : public std::runtime_error {

    public :
        explicit Exception( std::string const &a_message );
};

/*
============================================================
================== StatusMessageReporting ==================
============================================================
*/

class StatusMessageReporting {

    public:
        enum class Status { ok, info, warning, error };

    private:
        statusMessageReporting m_smr;

    public:
        StatusMessageReporting( );
        ~StatusMessageReporting( );
        
        statusMessageReporting *smr( ) { return( &m_smr ); }
        bool isOk( ) { return( smr_isOk( &m_smr ) ); }
        bool isInfo( ) { return( smr_isInfo( &m_smr ) ); }
        bool isWarning( ) { return( smr_isWarning( &m_smr ) ); }
        bool isError( ) { return( smr_isError( &m_smr ) ); }
        void clear( ) { smr_release( &m_smr ); }
        std::string constructMessage( std::string a_prefix, int a_reports = 1, bool a_clear = false );
        std::string constructFullMessage( std::string const &a_prefix, int a_reports = 1, bool a_clear = false );
};

/*
============================================================
====================== ArgumentParser ======================
============================================================
*/

enum class ArgumentType { True, False, Count, Store, Append, Positional };

class ArgumentBase;

class ArgumentParser {

    private:
        std::string m_codeName;                                 /**< The name of the code that is using **ArgumentParser**. */
        std::string m_descriptor;                               /**< The descriptor that is printed when help (i.e., '-h') is entered. */
        std::vector<ArgumentBase *> m_arguments;                /**< The list of arguments (positional and optional) supported. */

        void add2( ArgumentBase *a_argumentBase );

    public:
        ArgumentParser( std::string const &a_codeName, std::string const &a_descriptor = "" );
        ~ArgumentParser( );

        std::string const &codeName( ) const { return( m_codeName ); }
        std::string const &descriptor( ) const { return( m_descriptor ); }
        template<typename T> T *add( std::string const &a_name, std::string const &a_descriptor, int a_minimumNeeded = 1, int a_maximumNeeded = 1 );
        ArgumentBase *add( ArgumentType a_argumentType, std::string const &a_name, std::string const &a_descriptor, 
                int a_minimumNeeded = -2, int a_maximumNeeded = -2 );
        void addAlias( std::string const &a_name, std::string const &a_alias );
        void addAlias( ArgumentBase const * const a_argumentBase, std::string const &a_alias );
        bool hasName( std::string const &a_name ) const ;
        bool isOptionalArgument( std::string const &a_name ) const ;
        void parse( int a_argc, char **a_argv, bool a_printArguments = true );
        template<typename T> T *get( std::size_t a_name );
        void help( ) const ;
        void usage( ) const ;
        virtual void printStatus( std::string a_indent ) const ;
};

/* *********************************************************************************************************//**
 * Creates a new argument, adds the argument to *this* and returns a pointer the the newly created argument.
 *
 * @param a_name                [in]    The name of the argument.
 * @param a_descriptor          [in]    The argument's description, displayed when the help option is enetered.
 * @param a_minimumNeeded       [in]    The minimum number of required time *this* argument must be entered.
 * @param a_maximumNeeded       [in]    The maximum number of required time *this* argument must be entered.
 *
 * @return                              A pointer to the created argument.
 ***********************************************************************************************************/

template<typename T> T *ArgumentParser::add( std::string const &a_name, std::string const &a_descriptor, int a_minimumNeeded, int a_maximumNeeded ) {

    T *argument = new T( a_name, a_descriptor, a_minimumNeeded, a_maximumNeeded );
    add2( argument );

    return( argument );
}

/*
============================================================
======================= ArgumentBase =======================
============================================================
*/

class ArgumentBase {

    private:
        ArgumentType m_argumentType;                        /**< The enum for arguent type of *this*. */
        std::vector<std::string> m_names;                   /**< The allowed names for *this*. */
        std::string m_descriptor;                           /**< The desciption printed help. */
        int m_minimumNeeded;                                /**< Minimum number of times *this* argument is required on the command line. */
        int m_maximumNeeded;                                /**< Maximum number of times *this* argument is required on the command line. */
        std::size_t m_counts;                               /**< The number of time this argument was entered on the command line. */
        std::vector<std::string> m_values;                  /**< list of values entered for this argument. Only used for types Store, Append and Positional. */

        void addAlias( std::string const &a_name );                         /**< Adds the alias *a_name* to *this*. */
        virtual std::string printStatus2( ) const ;                         /**< For internal use. Called by method **printStatus**. */
        virtual void printStatus3( std::string const &a_indent ) const ;

        friend void ArgumentParser::addAlias( std::string const &a_name, std::string const &a_alias );

    public:
        ArgumentBase( ArgumentType a_argumentType, std::string const &a_name, std::string const &a_descriptor, int a_minimumNeeded, int a_maximumNeeded );
        virtual ~ArgumentBase( ) = 0 ;

        ArgumentType argumentType( ) const { return( m_argumentType ); }
        std::string const &name( ) const { return( m_names[0] ); }
        std::vector<std::string> const &names( ) { return( m_names ); }
        bool hasName( std::string const &a_name ) const ;
        std::string const &descriptor( ) const { return( m_descriptor ); }
        int minimumNeeded( ) const { return( m_minimumNeeded ); }
        int maximumNeeded( ) const { return( m_maximumNeeded ); }
        std::size_t counts( ) const { return( m_counts ); }

        virtual std::string const &value( std::size_t a_index = 0 ) const ;
        std::vector<std::string> const &values( ) const { return( m_values ); }
        virtual bool isOptionalArgument( ) const { return( true ); }
        virtual bool requiresAValue( ) const { return( false ); }
        virtual int parse( ArgumentParser const &a_argumentParser, int a_index, int a_argc, char **a_argv );
        std::string usage( bool a_requiredOption ) const ;
        void printStatus( std::string a_indent ) const ;
};

/*
============================================================
======================= OptionBoolean ======================
============================================================
*/

class OptionBoolean : public ArgumentBase {

    private:
        bool m_default;

    public:
        OptionBoolean( ArgumentType a_argumentType, std::string const &a_name, std::string const &a_descriptor, bool a_default );
        virtual ~OptionBoolean( ) = 0 ;

        bool _default( ) const { return( m_default ); }
        std::string printStatus2( ) const ;
};

/*
============================================================
======================== OptionTrue ========================
============================================================
*/

class OptionTrue : public OptionBoolean {

    public:
        OptionTrue( std::string const &a_name, std::string const &a_descriptor = "", int a_minimumNeeded = 0, int a_maximumNeeded = -1 );
        ~OptionTrue( ) { }
};

/*
============================================================
======================= OptionFalse ========================
============================================================
*/

class OptionFalse : public OptionBoolean {

    public:
        OptionFalse( std::string const &a_name, std::string const &a_descriptor = "", int a_minimumNeeded = 0, int a_maximumNeeded = -1 );
        ~OptionFalse( ) { }
};

/*
============================================================
====================== OptionCounter =======================
============================================================
*/

class OptionCounter : public ArgumentBase {

    public:
        OptionCounter( std::string const &a_name, std::string const &a_descriptor = "", int a_minimumNeeded = 0, int a_maximumNeeded = -1 );
        ~OptionCounter( ) { }

        std::string printStatus2( ) const ;
};

/*
============================================================
======================= OptionStore ========================
============================================================
*/

class OptionStore : public ArgumentBase {

    public:
        OptionStore( std::string const &a_name, std::string const &a_descriptor = "", int a_minimumNeeded = 0, int a_maximumNeeded = -1 );
        ~OptionStore( ) { }

        std::string const &value( std::size_t a_index = 0 ) const ;
        bool requiresAValue( ) const { return( true ); }
        void printStatus3( std::string const &a_indent ) const ;
};

/*
============================================================
======================= OptionAppend =======================
============================================================
*/

class OptionAppend : public ArgumentBase {

    public:
        OptionAppend( std::string const &a_name, std::string const &a_descriptor = "", int a_minimumNeeded = 0, int a_maximumNeeded = -1 );
        ~OptionAppend( ) { }

        bool requiresAValue( ) const { return( true ); }
        void printStatus3( std::string const &a_indent ) const ;
};

/*
============================================================
======================== Positional ========================
============================================================
*/

class Positional : public ArgumentBase {

    public:
        Positional( std::string const &a_name, std::string const &a_descriptor = "", int a_minimumNeeded = 1, int a_maximumNeeded = 1 );
        ~Positional( ) { }

        bool isOptionalArgument( ) const { return( false ); }
        bool requiresAValue( ) const { return( true ); }
        void printStatus3( std::string const &a_indent ) const ;
};

/*
============================================================
======================== DeltaTime =========================
============================================================
*/

#ifndef _WIN32

#define LUPI_DeltaTime_toStringFormatIncremental "incremental: CPU %8.3fs, wall %8.3fs"
#define LUPI_DeltaTime_toStringFormatTotal "total: CPU %8.3fs, wall %8.3fs"

class DeltaTime {

    private:
        double m_CPU_time;
        double m_wallTime;
        double m_CPU_timeIncremental;
        double m_wallTimeIncremental;

    public:
        DeltaTime( );
        DeltaTime( double a_CPU_time, double a_wallTime, double a_CPU_timeIncremental, double a_wallTimeIncremental );
        DeltaTime( DeltaTime const &deltaTime );
        ~DeltaTime( ) { }

        double CPU_time( ) const { return( m_CPU_time ); }
        double wallTime( ) const { return( m_wallTime ); }
        double CPU_timeIncremental( ) const { return( m_CPU_timeIncremental ); }
        double wallTimeIncremental( ) const { return( m_wallTimeIncremental ); }
        std::string toString( std::string a_formatIncremental = LUPI_DeltaTime_toStringFormatIncremental,
                std::string a_format = LUPI_DeltaTime_toStringFormatTotal, std::string a_sep = "; " );
};

/*
============================================================
========================== Timer ===========================
============================================================
*/

class Timer {

    private:
        clock_t m_CPU_time;
        struct timeval m_wallTime;
        clock_t m_CPU_timeIncremental;
        struct timeval m_wallTimeIncremental;

    public:
        Timer( );
        ~Timer( ) { }

        DeltaTime deltaTime( );
        DeltaTime deltaTimeAndReset( );
        void reset( );
};

#endif          // End of not _WIN32 defined.

namespace FileInfo {        // Should be using std::filesystem stuff but this requires C++ 17.

std::string realPath( std::string const &a_path );
std::string _basename( std::string const &a_path );
std::string basenameWithoutExtension( std::string const &a_path );
std::string _dirname( std::string const &a_path );
bool exists( std::string const &a_path );
bool isDirectory( std::string const &a_path );
bool createDirectories( std::string const &a_path );

/*
============================================================
========================= FileStat =========================
============================================================
*/
class FileStat {

    private:
        std::string m_path;             /**< The path that is stat-ed. */
        struct stat m_stat;             /**< The stat for the path. */

    public:
        FileStat( std::string const &a_path );

        std::string const &path( ) const { return( m_path ); }                              /**< Returns a reference to the **m_path** member. */
        struct stat const &statRef( ) const { return( m_stat ); }                           /**< Returns a reference to the **m_stat** member. */

        bool exists( );
        bool isDirectory( ) const { return( ( m_stat.st_mode & S_IFMT ) == S_IFDIR ); }           /**< Returns *true* if the path is a directory and *false* otherwise. */
        bool isRegularFile( ) const { return( ( m_stat.st_mode & S_IFMT ) == S_IFREG ); }   /**< Returns *true* if the path is a regular file and *false* otherwise. */
};

}               // End of namespace FileInfo.

// Miscellaneous functions

namespace Misc {

std::string stripString( std::string const &a_string, bool a_left = true, bool a_right = true );
std::vector<std::string> splitString( std::string const &a_string, char a_delimiter, bool a_strip = false );
std::vector<std::string> splitString( std::string const &a_string, std::string const &a_delimiter, bool a_strip = false );
std::string joinStrings( std::string const &a_sep, std::vector<std::string> a_strings );
std::string replaceString( std::string const &a_string, std::string const &a_old, std::string const &a_new, bool a_all );
std::vector<std::string> splitXLinkString( std::string const &a_string );
bool stringToInt( std::string const &a_string, int &a_value );
bool stringToSize_t( std::string const &a_string, std::size_t &a_value );

std::string argumentsToString( char const *a_format, ... );
std::string doubleToString3( char const *a_format, double a_value, bool a_reduceBits = false );
std::string doubleToShortestString( double a_value, int a_significantDigits = 15, int a_favorEFormBy = 0 );

void printCommand( std::string const &a_indent, int a_argc, char **a_argv );

}               // End of namespace Misc.

}               // End of namespace LUPI.

#endif          // LUPI_hpp_included
