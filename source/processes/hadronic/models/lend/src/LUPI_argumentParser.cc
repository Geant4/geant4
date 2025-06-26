/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <iostream>

#include <LUPI.hpp>

namespace LUPI {

static void printArgumentDescription( std::string const &a_line, std::string const &a_descriptor );
static std::size_t maxPrintLineWidth = 120;

/*! \class ArgumentBase
 * Base class for argument and option sub-classes.
 */

/* *********************************************************************************************************//**
 * ArgumentBase constructor.
 *
 * @param a_argumentType        [in]    The type of argument to create.
 * @param a_name                [in]    The name of the argument.
 * @param a_descriptor          [in]    The string printed with arugment's help.
 * @param a_minimumNeeded       [in]    The minimum number of times the argument must be entered.
 * @param a_maximumNeeded       [in]    The maximum number of times the argument can be entered.
 ***********************************************************************************************************/

ArgumentBase::ArgumentBase( ArgumentType a_argumentType, std::string const &a_name, std::string const &a_descriptor, int a_minimumNeeded, int a_maximumNeeded ) :
        m_argumentType( a_argumentType ),
        m_names( ),
        m_descriptor( a_descriptor ),
        m_minimumNeeded( a_minimumNeeded ),
        m_maximumNeeded( a_maximumNeeded ),
        m_counts( 0 ) {

    if( m_minimumNeeded < 0 ) throw std::runtime_error( "ERROR 1000 in ArgumentBase::ArgumentBase: m_minimumNeeded must not be negative." );
    if( m_maximumNeeded > -1 ) {
        if( m_minimumNeeded > m_maximumNeeded )
            throw std::runtime_error( "ERROR 1010 in ArgumentBase::ArgumentBase: for argument '" + a_name + "' m_maximumNeeded less than m_minimumNeeded." );
    }

    if( m_argumentType == ArgumentType::Positional ) {
        if( a_name[0] == '-' )
            throw std::runtime_error( "ERROR 1020 in ArgumentBase::ArgumentBase: positional argument name '" + a_name + "' cannot start with a '-'." );
    }

    addAlias( a_name );
}

/* *********************************************************************************************************//**
 * ArgumentBase destructor.
 ***********************************************************************************************************/

ArgumentBase::~ArgumentBase( ) {

}

/* *********************************************************************************************************//**
 * Returns true if *a_name* is one of the names for *this* and false otherwise.
 *
 * @param a_name            [in]    The name to search for.
 *
 * @return                          Returns true if *a_name* is a match for one of the names for *this* and false otherwise.
 ***********************************************************************************************************/

bool ArgumentBase::hasName( std::string const &a_name ) const {

    for( auto iter = m_names.begin( ); iter != m_names.end( ); ++iter ) {
        if( a_name == *iter ) return( true );
    }
    return( false );
}

/* *********************************************************************************************************//**
 * Returns the value at index *a_index*. If *a_index* exceeds the number of values entered, a throw is executed.
 *
 * @param a_index               [in]    The 0-based index into the m_values std::vector whose content is returned.
 *
 * @return                              Returns the value entered at index *a_index*.
 ***********************************************************************************************************/

std::string const &ArgumentBase::value( std::size_t a_index ) const {

    if( ( m_argumentType == ArgumentType::True ) || ( m_argumentType == ArgumentType::False ) || ( m_argumentType == ArgumentType::Count ) )
        throw Exception( "Argument type for " + name( ) + " does not support calling value() method." );

    if( a_index >= m_values.size( ) ) throw Exception( "Index = " + std::to_string( a_index ) + " out-of-bounds for argument \"" + name( ) + "\"." );

    return( m_values[a_index] );
}

/* *********************************************************************************************************//**
 * Add *a_name* as an optional name for *this*.
 *
 * @param a_name            [in]    The name to add.
 ***********************************************************************************************************/

void ArgumentBase::addAlias( std::string const &a_name ) {

    if( m_argumentType == ArgumentType::Positional ) {
        if( m_names.size( ) > 0 ) throw std::runtime_error( "ERROR 1100 in ArgumentBase::addAlias: cannot add a name to a positional argument." ); }
    else {
        if( a_name[0] != '-' ) throw std::runtime_error( "ERROR 1110 in ArgumentBase::addAlias: name '" + a_name + "' not a valid optional name." );
    }

    if( hasName( a_name ) ) return;

    m_names.push_back( a_name );
}

/* *********************************************************************************************************//**
 * Counts each time a specific argument is found.
 *
 * @param a_index           [in]    The index of the current command argument in *a_argv*.
 * @param a_argc            [in]    The number of command arguments.
 * @param a_argv            [in]    The list of command arguments.
 *
 * @return                          The value of *a_index* + 1.
 ***********************************************************************************************************/

int ArgumentBase::parse( ArgumentParser const &a_argumentParser, int a_index, int a_argc, char **a_argv ) {

    ++m_counts;
    if( m_argumentType != ArgumentType::Positional ) ++a_index;

    if( ( m_argumentType == ArgumentType::Store ) || ( m_argumentType == ArgumentType::Append ) || ( m_argumentType == ArgumentType::Positional ) ) {
        int maximumNeeded1 = maximumNeeded( );
        if( m_argumentType == ArgumentType::Positional ) {
            if( maximumNeeded1 < 0 ) maximumNeeded1 = a_argc; }
        else {
            if( ( maximumNeeded1 < counts( ) ) && ( maximumNeeded1 > -1 ) )
                throw std::runtime_error( "ERROR 1220 in ArgumentBase::parse: too many values for optional argument " + name( ) + " entered." );
            maximumNeeded1 = 1;
        }

        for( int index = 0; index < maximumNeeded1; ++index ) {
            if( a_index == a_argc ) {
                if( m_argumentType == ArgumentType::Positional ) break;
                throw std::runtime_error( "ERROR 1200 in ArgumentBase::parse: missing value for argument " + name( ) + "." );
            }
            if( ( m_argumentType == ArgumentType::Positional ) && a_argumentParser.isOptionalArgument( a_argv[a_index] ) ) break;

            m_values.push_back( a_argv[a_index] );
            if( index > 0 ) ++m_counts;
            ++a_index;
        }
    }

    return( a_index );
}

/* *********************************************************************************************************//**
 * Returns the usage string for *this* argument.
 *
 * @param a_requiredOption      [in]    The index of the current command argument in *a_argv*.
 *
 * @return                              The value of *a_index* + 1.
 ***********************************************************************************************************/

std::string ArgumentBase::usage( bool a_requiredOption ) const {

    std::string usageString;

    if( isOptionalArgument( ) ) {
        if( a_requiredOption ) {
            if( m_minimumNeeded != 0 ) return( usageString ); }
        else {
            if( m_minimumNeeded == 0 ) return( usageString );
        }
    }
    usageString += " ";

    std::string value;
    if( isOptionalArgument( ) && requiresAValue( ) ) value = " VALUE";
    if( isOptionalArgument( ) && ( m_minimumNeeded == 0 ) ) {
        usageString += "[" + name( ) + value + "]"; }
    else {
        usageString += name( ) + value;
    }

    if( !isOptionalArgument( ) ) {
        if( ( m_minimumNeeded != 1 ) || ( m_maximumNeeded != 1 ) )
            usageString += "[" + std::to_string( m_minimumNeeded ) + "," + std::to_string( m_maximumNeeded ) + "]";
    }

    return( usageString );
}

/* *********************************************************************************************************//**
 * Prints generic information about the status of *this*.
 *
 * @param a_indent          [in]    The amount of indentation to start the first line with.
 *
 * @return                          The value of *a_index* + 1.
 ***********************************************************************************************************/

void ArgumentBase::printStatus( std::string a_indent ) const {

    std::string name1 = name( );
    if( name1.size( ) < 32 ) name1.resize( 32, ' ' );

    std::cout << a_indent << name1 << ": number entered " << std::to_string( m_counts )
            << "; number needed (" << std::to_string( m_minimumNeeded ) << "," << std::to_string( m_maximumNeeded ) << ")" 
            << printStatus2( ) << std::endl;
    printStatus3( a_indent + "    " );
}

/* *********************************************************************************************************//**
 * Called by *printStatus*. This method returns an empty string. Must be overwritten by argument classes that have value.
 *
 * @return                          An empty **std::string** instance.
 ***********************************************************************************************************/

std::string ArgumentBase::printStatus2( ) const {

    return( "" );

}

/* *********************************************************************************************************//**
 * Called by *printStatus*. This method does nothing. Must be overwritten by argument classes that have value(s).
 *
 * @param a_indent          [in]    The amount of indentation to start the first line with.
 ***********************************************************************************************************/

void ArgumentBase::printStatus3( LUPI_maybeUnused std::string const &a_indent ) const {

}

/*! \class OptionBoolean
 * Base boolean class.
 */

/* *********************************************************************************************************//**
 * OptionBoolean constructor.
 *
 * @param a_argumentType        [in]    The type of argument to create.
 * @param a_name                [in]    The name of the argument.
 * @param a_descriptor          [in]    The string printed with arugment's help.
 * @param a_default             [in]    The default bool value.
 ***********************************************************************************************************/

OptionBoolean::OptionBoolean( ArgumentType a_argumentType, std::string const &a_name, std::string const &a_descriptor, bool a_default ) :
        ArgumentBase( a_argumentType, a_name, a_descriptor, 0, -1 ),
        m_default( a_default ) {

}

/* *********************************************************************************************************//**
 * OptionBoolean destructor.
 ***********************************************************************************************************/

OptionBoolean::~OptionBoolean( ) {

}

/* *********************************************************************************************************//**
 * Called by *printStatus*. This method returns a string representing *this*'s value.
 *
 * @return                              Returns a std::string instance representing the value of *this*.
 ***********************************************************************************************************/

std::string OptionBoolean::printStatus2( ) const {

    bool value1 = m_default;
    if( counts() == 0 ) value1 = !m_default;

    if( value1 ) return( ": true" );
    return( ": false" );
}

/*! \class OptionTrue
 * An boolean optional argument whose default is *false* and changes to *true* if one or more options are entered.
 */

/* *********************************************************************************************************//**
 * OptionTrue constructor.
 *
 * @param a_name                [in]    The name of the argument.
 * @param a_descriptor          [in]    The string printed with arugment's help.
 * @param a_minimumNeeded       [in]    Not used. Will probably be deprecated.
 * @param a_maximumNeeded       [in]    Not used. Will probably be deprecated.
 ***********************************************************************************************************/

OptionTrue::OptionTrue( std::string const &a_name, std::string const &a_descriptor, LUPI_maybeUnused int a_minimumNeeded, LUPI_maybeUnused int a_maximumNeeded ) :
        OptionBoolean( ArgumentType::True, a_name, a_descriptor, false ) {

}

/*! \class OptionFalse
 * An boolean optional argument whose default is *true* and changes to *false* if one or more options are entered.
 */

/* *********************************************************************************************************//**
 * OptionFalse constructor.
 *
 * @param a_name                [in]    The name of the argument.
 * @param a_descriptor          [in]    The string printed with arugment's help.
 * @param a_minimumNeeded       [in]    Not used. Will probably be deprecated.
 * @param a_maximumNeeded       [in]    Not used. Will probably be deprecated.
 ***********************************************************************************************************/

OptionFalse::OptionFalse( std::string const &a_name, std::string const &a_descriptor, LUPI_maybeUnused int a_minimumNeeded, LUPI_maybeUnused int a_maximumNeeded ) :
        OptionBoolean( ArgumentType::False, a_name, a_descriptor, true ) {

}

/*! \class OptionCounter
 * An optional argument that counts the number of times the option is entered.
 */

/* *********************************************************************************************************//**
 * OptionCounter constructor.
 *
 * @param a_name                [in]    The name of the argument.
 * @param a_descriptor          [in]    The string printed with arugment's help.
 * @param a_minimumNeeded       [in]    Not used. Will probably be deprecated.
 * @param a_maximumNeeded       [in]    Not used. Will probably be deprecated.
 ***********************************************************************************************************/

OptionCounter::OptionCounter( std::string const &a_name, std::string const &a_descriptor, LUPI_maybeUnused int a_minimumNeeded, LUPI_maybeUnused int a_maximumNeeded ) :
        ArgumentBase( ArgumentType::Count, a_name, a_descriptor, 0, -1 ) {

}

/* *********************************************************************************************************//**
 * Called by *printStatus*. This method returns a string representing *this*'s value.
 ***********************************************************************************************************/

std::string OptionCounter::printStatus2( ) const {

    return( " counts = " + std::to_string( counts( ) ) );
}

/*! \class OptionStore
 * An option with a value. If multiple options with the same name are entered, only the last vluae entered if returned by the value() method.
 */

/* *********************************************************************************************************//**
 * OptionStore constructor.
 *
 * @param a_name                [in]    The name of the argument.
 * @param a_descriptor          [in]    The string printed with arugment's help.
 * @param a_minimumNeeded       [in]    Not used. Will probably be deprecated.
 * @param a_maximumNeeded       [in]    Not used. Will probably be deprecated.
 ***********************************************************************************************************/

OptionStore::OptionStore( std::string const &a_name, std::string const &a_descriptor, LUPI_maybeUnused int a_minimumNeeded, LUPI_maybeUnused int a_maximumNeeded ) :
        ArgumentBase( ArgumentType::Store, a_name, a_descriptor, 0, -1 ) {

}

/* *********************************************************************************************************//**
 * Returns the value of 
 *
 * @param a_index               [in]    This argument is not used. The last value entered is always returned.
 *
 * @return                              Returns the last value entered or executes a **throw** if option not entered.
 ***********************************************************************************************************/

std::string const &OptionStore::value( std::size_t a_index ) const {

    a_index = values( ).size( );
    if( a_index != 0 ) --a_index;

    return ArgumentBase::value( a_index );
}

/* *********************************************************************************************************//**
 * Prints the value for *this*. Called by *printStatus*. 
 *
 * @param a_indent          [in]    The amount of indentation to start the first line with.
 ***********************************************************************************************************/

void OptionStore::printStatus3( std::string const &a_indent ) const {

    if( counts( ) > 0 ) std::cout << a_indent << value( ) << std::endl;
}

/*! \class OptionAppend
 * An option with a value. If multiple options with the same name are entered, all values will be stored.
 */

/* *********************************************************************************************************//**
 * OptionAppend constructor.
 *
 * @param a_name                [in]    The name of the argument.
 * @param a_descriptor          [in]    The string printed with arugment's help.
 * @param a_minimumNeeded       [in]    The minimum number of times the argument must be entered.
 * @param a_maximumNeeded       [in]    The maximum number of times the argument can be entered.
 ***********************************************************************************************************/

OptionAppend::OptionAppend( std::string const &a_name, std::string const &a_descriptor, int a_minimumNeeded, int a_maximumNeeded ) :
        ArgumentBase( ArgumentType::Append, a_name, a_descriptor, a_minimumNeeded, a_maximumNeeded ) {

}

/* *********************************************************************************************************//**
 * Prints the values for *this*. Called by *printStatus*. 
 *
 * @param a_indent          [in]    The amount of indentation to start the first line with.
 ***********************************************************************************************************/

void OptionAppend::printStatus3( std::string const &a_indent ) const {

    for( auto valueIterator = values( ).begin( ); valueIterator != values( ).end( ); ++valueIterator ) {
        std::cout << a_indent << *valueIterator << std::endl;
    }
}

/*! \class Positional
 * An option with a value. If multiple options with the same name are entered, the *m_value* member will represent the last option entered.
 */ 

/* *********************************************************************************************************//**
 * Positional constructor.
 *
 * @param a_name                [in]    The name of the argument.
 * @param a_descriptor          [in]    The string printed with arugment's help.
 * @param a_minimumNeeded       [in]    The minimum number of times the argument must be entered.
 * @param a_maximumNeeded       [in]    The maximum number of times the argument can be entered.
 ***********************************************************************************************************/

Positional::Positional( std::string const &a_name, std::string const &a_descriptor, int a_minimumNeeded, int a_maximumNeeded ) :
            ArgumentBase( ArgumentType::Positional, a_name, a_descriptor, a_minimumNeeded, a_maximumNeeded ) {
}

/* *********************************************************************************************************//**
 * Prints the values for *this*. Called by *printStatus*. 
 *
 * @param a_indent          [in]    The amount of indentation to start the first line with.
 ***********************************************************************************************************/

void Positional::printStatus3( std::string const &a_indent ) const {

    for( auto valueIterator = values( ).begin( ); valueIterator != values( ).end( ); ++valueIterator ) {
        std::cout << a_indent << *valueIterator << std::endl;
    }
}

/*! \class ArgumentParser
 * The main argument parser class.
 */

/* *********************************************************************************************************//**
 * ArgumentParser constructor.
 ***********************************************************************************************************/

ArgumentParser::ArgumentParser( std::string const &a_codeName, std::string const &a_descriptor ) :
        m_codeName( FileInfo::basenameWithoutExtension( a_codeName ) ),
        m_descriptor( a_descriptor ) {

}

/* *********************************************************************************************************//**
 * ArgumentParser destructor.
 ***********************************************************************************************************/

ArgumentParser::~ArgumentParser( ) {

    for( auto argumentIterator = m_arguments.begin( ); argumentIterator != m_arguments.end( ); ++argumentIterator ) {
        delete *argumentIterator;
    }
}

/* *********************************************************************************************************//**
 * Add *a_argumentBase* to the list of arguments.
 *
 * @param a_argumentBase   [in]     Pointer to the *ArgumentBase* to add to *this*.
 ***********************************************************************************************************/

void ArgumentParser::add2( ArgumentBase *a_argumentBase ) {

    if( !a_argumentBase->isOptionalArgument( ) ) {
        for( auto argumentIterator = m_arguments.rbegin( ); argumentIterator != m_arguments.rend( ); ++argumentIterator ) {
            if( !(*argumentIterator)->isOptionalArgument( ) ) {
                if( (*argumentIterator)->minimumNeeded( ) == (*argumentIterator)->maximumNeeded( ) ) break;
                throw std::runtime_error( "ERROR 1400 in ArgumentParser::add: request to add postional argument when prior postional argument '" 
                        + (*argumentIterator)->name( ) + "' takes a variable number of values." );
            }
        }
    }

    if( hasName( a_argumentBase->name( ) ) )
        throw std::runtime_error( "ERROR 1500 in ArgumentParser::add: name '" + a_argumentBase->name( ) + "' already present." );

    m_arguments.push_back( a_argumentBase );
}

/* *********************************************************************************************************//**
 * Creates an argument instance of type specified by *a_argumentType* and adds to *this*.
 *
 * @param a_argumentType        [in]    The type of argument to create.
 * @param a_name                [in]    The name of the argument.
 * @param a_descriptor          [in]    The string printed with arugment's help.
 * @param a_minimumNeeded       [in]    The minimum number of times the argument must be entered.
 * @param a_maximumNeeded       [in]    The maximum number of times the argument can be entered.
 *
 * @return                              Returns a pointer to to the created argument instance.
 ***********************************************************************************************************/

ArgumentBase *ArgumentParser::add( ArgumentType a_argumentType, std::string const &a_name, std::string const &a_descriptor, 
                int a_minimumNeeded, int a_maximumNeeded ) {

    ArgumentBase *argument = nullptr;

    switch( a_argumentType ) {
    case ArgumentType::True :
        argument = new OptionTrue( a_name, a_descriptor );
        break;
    case ArgumentType::False :
        argument = new OptionFalse( a_name, a_descriptor );
        break;
    case ArgumentType::Count :
        argument = new OptionCounter( a_name, a_descriptor );
        break;
    case ArgumentType::Store :
        argument = new OptionStore( a_name, a_descriptor );
        break;
    case ArgumentType::Append :
        argument = new OptionAppend( a_name, a_descriptor, a_minimumNeeded, a_maximumNeeded );
        break;
    default :
        argument = new Positional( a_name, a_descriptor, a_minimumNeeded, a_maximumNeeded );
    }

    add2( argument );

    return( argument );
}

/* *********************************************************************************************************//**
 * Adds the alias *a_alias* to the argument named *a_name*.
 *
 * @param a_name            [in]    The name of the argument to add the alias to.
 * @param a_alias           [in]    The alias name to add.
 ***********************************************************************************************************/

void ArgumentParser::addAlias( std::string const &a_name, std::string const &a_alias ) {

    if( hasName( a_alias ) )
        throw std::runtime_error( "ERROR 1510 in ArgumentParser::addAlias: name '" + a_alias + "' already present." );

    for( auto argumentIterator = m_arguments.begin( ); argumentIterator != m_arguments.end( ); ++argumentIterator ) {
        if( (*argumentIterator)->hasName( a_name ) ) {
            (*argumentIterator)->addAlias( a_alias );
            return;
        }
    }
    throw std::runtime_error( "ERROR 1520 in ArgumentParser::addAlias: no such argument named '" + a_name + "'." );
}

/* *********************************************************************************************************//**
 * Adds the alias *a_alias* to the argument *a_argumentBase*.
 *
 * @param a_argumentBase    [in]    The argument to add the alias to.
 * @param a_alias           [in]    The name of the argument to add the alias to.
 ***********************************************************************************************************/

void ArgumentParser::addAlias( ArgumentBase const * const a_argumentBase, std::string const &a_alias ) {

    addAlias( a_argumentBase->name( ), a_alias );
}

/* *********************************************************************************************************//**
 * Returns true if name *a_name* is an optional argument of *this* and false otherwise.
 *
 * @param a_name            [in]    The name to see check if it exists in *this*.
 *
 * @return                          true if name *a_name* is in *this* and false otherwise.
 ***********************************************************************************************************/

bool ArgumentParser::isOptionalArgument( std::string const &a_name ) const {

    for( auto argumentIterator = m_arguments.begin( ); argumentIterator != m_arguments.end( ); ++argumentIterator ) {
        if( (*argumentIterator)->hasName( a_name ) ) return( (*argumentIterator)->argumentType( ) != ArgumentType::Positional );
    }

    return( false );
}

/* *********************************************************************************************************//**
 * Returns true if name *a_name* is in *this* and false otherwise.
 *
 * @param a_name            [in]    The name to see check if it exists in *this*.
 *
 * @return                          true if name *a_name* is in *this* and false otherwise.
 ***********************************************************************************************************/

bool ArgumentParser::hasName( std::string const &a_name ) const {

    for( auto argumentIterator = m_arguments.begin( ); argumentIterator != m_arguments.end( ); ++argumentIterator ) {
        if( (*argumentIterator)->hasName( a_name ) ) return( true );
    }

    return( false );
}

/* *********************************************************************************************************//**
 * Parses the list of arguments.
 *
 * @param a_argc            [in]    The number of arguments.
 * @param a_argv            [in]    The list of arguments.
 ***********************************************************************************************************/

void ArgumentParser::parse( int a_argc, char **a_argv, bool a_printArguments ) {

    for( int iargc = 1; iargc < a_argc; ++iargc ) {                                       // Check is help requested.
        std::string arg( a_argv[iargc] );

        if( ( arg == "-h" ) || ( arg == "--help" ) ) help( );
    }

    auto argumentIterator = m_arguments.begin( );      // Find first non-option argument.
    for( ; argumentIterator != m_arguments.end( ); ++argumentIterator ) {
        if( !(*argumentIterator)->isOptionalArgument( ) ) break;
    }

    int iargc = 1;
    for( ; iargc < a_argc; ) {
        std::string arg( a_argv[iargc] );

        if( arg[0] == '-' ) {                                                           // Need to check if negative number.
            auto argumentIterator2 = m_arguments.begin( );
            for( ; argumentIterator2 != m_arguments.end( ); ++argumentIterator2 ) {
                if( (*argumentIterator2)->hasName( arg ) ) break;
            }
            if( argumentIterator2 == m_arguments.end( ) ) throw std::runtime_error( "ERROR 1600 in ArgumentParser::parse: invalid option '" + arg + "'." );
            iargc = (*argumentIterator2)->parse( *this, iargc, a_argc, a_argv ); }
        else {
            if( argumentIterator == m_arguments.end( ) )
                throw std::runtime_error( "ERROR 1610 in ArgumentParser::parse: additional positional argument found starting at index " 
                        + std::to_string( iargc ) + " (" + arg + ")" );

            iargc = (*argumentIterator)->parse( *this, iargc, a_argc, a_argv );

            ++argumentIterator;
            for( ; argumentIterator != m_arguments.end( ); ++argumentIterator ) {       // Find next positional arguments.
                if( !(*argumentIterator)->isOptionalArgument( ) ) break;
            }
        }
    }
    for( auto argumentIterator2 = m_arguments.begin( ); argumentIterator2 != m_arguments.end( ); ++argumentIterator2 ) {
        if( (*argumentIterator2)->counts( ) < (*argumentIterator2)->minimumNeeded( ) ) {
            std::string msg( "arguments for" );

            if( (*argumentIterator2)->isOptionalArgument( ) ) msg = "number of option";
            throw std::runtime_error( "ERROR 1620 in ArgumentParser::parse: insufficient " + msg + " '" + (*argumentIterator2)->name( ) 
                    + "' entered. Range of " + std::to_string( (*argumentIterator2)->minimumNeeded( ) ) + " to " 
                    + std::to_string( (*argumentIterator2)->maximumNeeded( ) ) + " required, " 
                    + std::to_string( (*argumentIterator2)->counts( ) ) + " entered." );
        }
    }

    if( a_printArguments ) {
        std::cerr << "    " << LUPI::FileInfo::basenameWithoutExtension( m_codeName );
        for( int i1 = 1; i1 < a_argc; i1++ ) std::cerr << " " << a_argv[i1];
        std::cerr << std::endl;
    }
}

/* *********************************************************************************************************//**
 * Prints the help for *this*.
 ***********************************************************************************************************/

void ArgumentParser::help( ) const {

    usage( );

    if( m_descriptor != "" ) {
        std::cout << std::endl << "Description:" << std::endl;
        std::cout << "    " << m_descriptor << std::endl;
    }

    bool printHeader = true;
    for( auto argumentIterator = m_arguments.begin( ); argumentIterator != m_arguments.end( ); ++argumentIterator ) {
        if( (*argumentIterator)->isOptionalArgument( ) ) continue;

        if( printHeader ) std::cout << std::endl << "positional arguments:" << std::endl;
        printHeader = false;

        std::string line = (*argumentIterator)->name( );
        if( ( (*argumentIterator)->minimumNeeded( ) != (*argumentIterator)->maximumNeeded( ) ) || ( (*argumentIterator)->maximumNeeded( ) != 1 ) )
                    line += " [" + std::to_string( (*argumentIterator)->minimumNeeded( ) ) + "," + std::to_string( (*argumentIterator)->maximumNeeded( ) ) + "]";
        printArgumentDescription( line, (*argumentIterator)->descriptor( ) );
    }

    std::cout << std::endl << "optional arguments:" << std::endl;
    std::cout << "  -h, --help                  Show this help message and exit." << std::endl;
    for( auto argumentIterator = m_arguments.begin( ); argumentIterator != m_arguments.end( ); ++argumentIterator ) {
        if( !(*argumentIterator)->isOptionalArgument( ) ) continue;

        std::string line;
        std::string sep;
        for( auto namesIterator = (*argumentIterator)->names( ).begin( ); namesIterator != (*argumentIterator)->names( ).end( ); ++namesIterator ) {
            line += sep + *namesIterator;
            sep = ", ";
        }
        if( (*argumentIterator)->requiresAValue( ) ) line += " VALUE";
        if( (*argumentIterator)->argumentType( ) == ArgumentType::Append ) {
            if( ( (*argumentIterator)->minimumNeeded( ) != (*argumentIterator)->maximumNeeded( ) ) || ( (*argumentIterator)->maximumNeeded( ) != 1 ) )
                    line += " [" + std::to_string( (*argumentIterator)->minimumNeeded( ) ) + "," + std::to_string( (*argumentIterator)->maximumNeeded( ) ) + "]";
        }
        printArgumentDescription( line, (*argumentIterator)->descriptor( ) );
    }

    exit( EXIT_SUCCESS );
}

/* *********************************************************************************************************//**
 * Prints the usage for *this*.
 ***********************************************************************************************************/

void ArgumentParser::usage( ) const {

    std::string line( "usage: " );
    line += codeName( );
    std::string indent( "" );
    indent.resize( line.size( ), ' ' );

    for( int counter = 0; counter < 2; ++counter ) {
        for( auto argumentIterator = m_arguments.begin( ); argumentIterator != m_arguments.end( ); ++argumentIterator ) {
            if( (*argumentIterator)->isOptionalArgument( ) ) {
                std::string optionUsage( (*argumentIterator)->usage( counter == 0 ) );

                if( ( line.size( ) + optionUsage.size( ) ) > maxPrintLineWidth ) {
                    std::cout << line << std::endl;
                    line = indent;
                }
                line += optionUsage;
            }
        }
    }

    for( auto argumentIterator = m_arguments.begin( ); argumentIterator != m_arguments.end( ); ++argumentIterator ) {
        if( (*argumentIterator)->isOptionalArgument( ) ) continue;
        std::string optionUsage( (*argumentIterator)->usage( false ) );

        if( ( line.size( ) + optionUsage.size( ) ) > maxPrintLineWidth ) {
            std::cout << line << std::endl;
            line = indent;
        }
        line += optionUsage;
    }
    if( line.size( ) > indent.size( ) ) std::cout << line << std::endl;
    std::cout << std::endl;
}

/* *********************************************************************************************************//**
 * Returns the usage string for *this* option.
 *
 * @param a_indent          [in]    The amount of indentation to start the first line with.
 *
 * @return                          The value of *a_index* + 1.
 ***********************************************************************************************************/

void ArgumentParser::printStatus( std::string a_indent ) const {

    for( auto argumentIterator = m_arguments.begin( ); argumentIterator != m_arguments.end( ); ++argumentIterator ) {
        (*argumentIterator)->printStatus( a_indent );
    }
}


/* *********************************************************************************************************//**
 * For internal use only.
 *
 * @param a_line                [in]    A string containing the help line for an argument up to the description string.
 * @param a_descriptor          [in]    The help description string.
 ***********************************************************************************************************/

static void printArgumentDescription( std::string const &a_line, std::string const &a_descriptor ) {

    std::string newLineSpaces = "                            ";
    auto size = newLineSpaces.size( );
    auto maxDescriptionWidth = maxPrintLineWidth - size;
    std::string line = "  " + a_line;
    if( line.size( ) < size ) line.resize( size, ' ' );
    std::cout << line;
    if( line.size( ) > size ) std::cout << std::endl << newLineSpaces;

    std::string descriptor = a_descriptor;
    while( descriptor.size( ) > 0 ) {
        auto length = descriptor.size( );
        if( length > maxDescriptionWidth ) {
            length = maxDescriptionWidth;
            length = descriptor.rfind( ' ', length );
            if( length == 0 ) length = descriptor.find( ' ', length );
        }
    
        std::cout << "  " << descriptor.substr( 0, length ) << std::endl;
        descriptor = descriptor.substr( length );
        auto firstNotOf = descriptor.find_first_not_of( " " );
        if( firstNotOf != descriptor.npos ) descriptor = descriptor.substr( firstNotOf );
        if( descriptor.size( ) > 0 ) std::cout << newLineSpaces;
    }
}

}               // End of namespace LUPI.
