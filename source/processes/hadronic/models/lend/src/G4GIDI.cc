//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include <G4GIDI.hh>

PoPI::Database G4GIDI_pops;
static bool G4GIDI_pops_initialized = false;

/* *********************************************************************************************************//**
 * Returns the familiar name for the projectile for *a_ip*.
 *
 * @param a_ip              [in]    One of the following ids for the projectile (0:photon, 1:neutron, 2:proton, 3:deutron, 4:triton, 5:helion or 6:alpha).
 ***********************************************************************************************************/

static std::string projectileStringFromID( int a_ip ) {

    switch( a_ip ) {
    case 0:
        return( PoPI::IDs::photon );
    case 1:
        return( PoPI::IDs::neutron );
    case 2:
        return( PoPI::IDs::proton );
    case 3:
        return( PoPI::IDs::familiarDeuteron );
    case 4:
        return( PoPI::IDs::familiarTriton );
    case 5:
        return( PoPI::IDs::familiarHelion );
    case 6:
        return( PoPI::IDs::familiarAlpha );
    default:
        throw LUPI::Exception( "Invalid projectile id " + std::to_string( a_ip ) );
    }

    return( "" );
}

/* *********************************************************************************************************//**
 * Initialize the G4GIDI stuff.
 ***********************************************************************************************************/

void G4GIDI_initialize( std::string const &a_dataPath ) {

    if( G4GIDI_pops_initialized ) return;

    G4GIDI_pops_initialized = true;
    G4GIDI_pops.addFile( a_dataPath + "/" + "pops.xml", false );
    G4GIDI_pops.addFile( a_dataPath + "/" + "metastables_alias.xml", false );
}

/*! \class G4GIDI
 * A class to store map files for a particular projectile.
 */

/* *********************************************************************************************************//**
 * @param a_ip              [in]    One of the following ids for the projectile (0:photon, 1:neutron, 2:proton, 3:deutron, 4:triton, 5:helion or 6:alpha).
 * @param a_dataDirectory   [in]    A path to a map file to load.
 ***********************************************************************************************************/

G4GIDI::G4GIDI( int a_ip, std::string const &a_dataDirectory ) :
        m_projectileIP( a_ip ),
        m_projectile( projectileStringFromID( a_ip ) ) {

    addDataDirectory( a_dataDirectory );
}

/* *********************************************************************************************************//**
 *
 * @param a_id              [in]    This argument is ignored but needed for backwards compatibility.
 * @param a_dataDirectories [in]    A list of paths to a map files to load.
 ***********************************************************************************************************/

G4GIDI::G4GIDI( int a_ip, std::list<std::string> const &a_dataDirectories ) :
        m_projectileIP( a_ip ),
        m_projectile( projectileStringFromID( a_ip ) ) {

    for( auto mapIter = a_dataDirectories.begin( ); mapIter != a_dataDirectories.end( ); ++mapIter ) {
        addDataDirectory( *mapIter );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

G4GIDI::~G4GIDI( ) {

    for( auto mapIter = m_maps.begin( ); mapIter != m_maps.end( ); ++mapIter ) delete *mapIter;
    for( auto protareIter = m_protares.begin( ); protareIter != m_protares.end( ); ++protareIter ) delete *protareIter;
}

/* *********************************************************************************************************//**
 * Adds the map file *a_dataDirectory* to *this*.
 *
 * @param a_dataDirectory   [in]    A path to a map file to load.
 ***********************************************************************************************************/

int G4GIDI::addDataDirectory( std::string const &a_dataDirectory ) {

    for( auto mapIter = m_maps.begin( ); mapIter != m_maps.end( ); ++mapIter ) {
        if( (*mapIter)->fileName( ) == a_dataDirectory ) return( 0 );
    }

    m_maps.push_back( new GIDI::Map::Map( a_dataDirectory, G4GIDI_pops ) );

    return( 0 );
}

/* *********************************************************************************************************//**
 * Removes the map file with path *a_dataDirectory* from *this*.
 *
 * @param a_dataDirectory   [in]    A path to a map file to unload.
 ***********************************************************************************************************/

int G4GIDI::removeDataDirectory( std::string const &a_dataDirectory ) {

    std::size_t index = 0;

    for( auto mapIter = m_maps.begin( ); mapIter!= m_maps.end( ); ++mapIter, ++index ) {
        if( a_dataDirectory == (*mapIter)->fileName( ) ) {
            delete *mapIter;
   
            ++index;         
            for( ; index < m_maps.size( ); ++index )  m_maps[index-1] = m_maps[index];
            m_maps[index-1] = nullptr;
            m_maps.resize( m_maps.size( ) - 1 );
        }
    }

    return( 0 );
}

/* *********************************************************************************************************//**
 * Removes the map file with path *a_dataDirectory* from *this*.
 *
 * @param a_dataDirectory   [in]    A path to a map file to unload.
 ***********************************************************************************************************/

std::string const G4GIDI::getDataDirectoryAtIndex( int a_index ) const {

    std::string nullString;

    if( ( a_index < 0 ) || ( a_index >= static_cast<int>( m_maps.size( ) ) ) ) return( nullString );

    return( m_maps[a_index]->fileName( ) );
}

/* *********************************************************************************************************//**
 * Returns a list of paths to all loaded map files.
 *
 * @return          Returns a pointer that the user owns (i.e., the user is responsible for free-ing).
 ***********************************************************************************************************/

std::vector<std::string> *G4GIDI::getDataDirectories( ) const {

    std::vector<std::string> *list = new std::vector<std::string>( );

    for( auto mapIter = m_maps.begin( ); mapIter!= m_maps.end( ); ++mapIter ) {
        list->push_back( (*mapIter)->fileName( ) );
    }

    return( list );
}

/* *********************************************************************************************************//**
 * Returns **true** if the specified target exists in a map file and **false** otherwise.
 *
 * @param a_lib_name            [in]    The evaluation. Call be an empty string.
 * @param a_Z                   [in]    The atomic number of the target.
 * @param a_A                   [in]    The mass number of the taret.
 * @param a_M                   [in]    The meta-stable index of the target.
 *
 * @return                              Boolean value.
 ***********************************************************************************************************/

bool G4GIDI::isThisDataAvailable( std::string const &a_lib_name, int a_Z, int a_A, int a_M ) const {

    return( isThisDataAvailable( a_lib_name, G4GIDI_Misc_Z_A_m_ToName( a_Z, a_A, a_M ) ) );
}

/* *********************************************************************************************************//**
 * Returns **true** if the specified target exists in a map file and **false** otherwise.
 *
 * @param a_lib_name            [in]    The evaluation. Call be an empty string.
 * @param a_targetName          [in]    The target PoPs id.
 *
 * @return                              Boolean value.
 ***********************************************************************************************************/

bool G4GIDI::isThisDataAvailable( std::string const &a_lib_name, std::string const &a_targetName ) const {

    for( auto mapIter = m_maps.begin( ); mapIter!= m_maps.end( ); ++mapIter ) {
        if( (*mapIter)->isProtareAvailable( m_projectile, a_targetName, "", a_lib_name ) ) return( true );
    }

    return( false );
}

/* *********************************************************************************************************//**
 * Returns the file path to the specified target or **nullptr** if the target does not exists.
 *
 * @param a_lib_name            [in]    The evaluation. Call be an empty string.
 * @param a_Z                   [in]    The atomic number of the target.
 * @param a_A                   [in]    The mass number of the taret.
 * @param a_M                   [in]    The meta-stable index of the target.
 *
 * @return                              A std::string of the protare's path.
 ***********************************************************************************************************/

std::string G4GIDI::dataFilename( std::string const &a_lib_name, int a_Z, int a_A, int a_M ) const {

    return( dataFilename( a_lib_name, G4GIDI_Misc_Z_A_m_ToName( a_Z, a_A, a_M ) ) );
}

/* *********************************************************************************************************//**
 * Returns the file path to the specified target or **nullptr** if the target does not exists.
 *
 * @param a_lib_name            [in]    The evaluation. Call be an empty string.
 * @param a_targetName          [in]    The target PoPs id.
 *
 * @return                              Boolean value.
 ***********************************************************************************************************/

std::string G4GIDI::dataFilename( std::string const &a_lib_name, std::string const &a_targetName ) const {

    for( auto mapIter = m_maps.begin( ); mapIter!= m_maps.end( ); ++mapIter ) {
        std::string path = (*mapIter)->protareFilename( m_projectile, a_targetName, "", a_lib_name );
        if( path != "" ) return( path );
    }

    return( "" );
}

/* *********************************************************************************************************//**
 * Determines target name from *a_Z*, *a_A* and *a_M* and calls **getNamesOfAvailableLibraries** with target's name, and returns
 * its return value.
 *
 * @param a_Z                   [in]    The atomic number of the target.
 * @param a_A                   [in]    The mass number of the taret.
 * @param a_M                   [in]    The meta-stable index of the target.
 *
 * @return                              Pointer to vector<string>.
 ***********************************************************************************************************/

std::vector<std::string> *G4GIDI::getNamesOfAvailableLibraries( G4int a_Z, G4int a_A, G4int a_M ) const {

    return( getNamesOfAvailableLibraries( G4GIDI_Misc_Z_A_m_ToName( a_Z, a_A, a_M ) ) );
}

/* *********************************************************************************************************//**
 * Returns the list of all evaluations that have a target named *a_targetName*. User is responsible for freeing the returned list.
 *
 * @return                              Pointer to vector<string>.
 ***********************************************************************************************************/

std::vector<std::string> *G4GIDI::getNamesOfAvailableLibraries( std::string const &a_targetName ) const {

    std::vector<std::string> *listOfLibraries = new std::vector<std::string>( );

    for( auto mapIter = m_maps.cbegin( ); mapIter != m_maps.cend( ); ++mapIter ) {
        std::vector<GIDI::Map::ProtareBase const *> entries;
        (*mapIter)->findProtareEntries( entries, std::regex( m_projectile ), std::regex( a_targetName ) );
        for( auto entryIter = entries.begin( ); entryIter != entries.end( ); ++entryIter ) {
            listOfLibraries->push_back( (*entryIter)->evaluation( ) );
        }
    }

    return( listOfLibraries );
}

/* *********************************************************************************************************//**
 * Returns the list of available targets for *this*.
 *
 * @return          Returns a pointer that the user owns (i.e., the user is responsible for free-ing).
 ***********************************************************************************************************/

std::vector<std::string> *G4GIDI::getNamesOfAvailableTargets( ) const {

    std::vector<std::string> *list = new std::vector<std::string>( );
    std::set<std::string> targetIDs;

    for( auto mapIter = m_maps.begin( ); mapIter!= m_maps.end( ); ++mapIter ) {
        auto protareBases = (*mapIter)->directory( m_projectile );

        for( auto protareBaseIter = protareBases.begin( ); protareBaseIter != protareBases.end( ); ++protareBaseIter ) {
            targetIDs.insert( (*protareBaseIter)->targetID( ) );
        }
    }

    for( auto targetIter = targetIDs.begin( ); targetIter != targetIDs.end( ); ++targetIter ) {
        (*list).push_back( (*targetIter) );
    }

    return( list );
}

/* *********************************************************************************************************//**
 * Returns the specified target or **nullptr** if the target does not exists.
 *
 * @param a_lib_name            [in]    The evaluation. Call be an empty string.
 * @param a_Z                   [in]    The atomic number of the target.
 * @param a_A                   [in]    The mass number of the taret.
 * @param a_M                   [in]    The meta-stable index of the target.
 *
 * @return                              Boolean value.
 ***********************************************************************************************************/

G4GIDI_target *G4GIDI::readTarget( std::string const &a_lib_name, int a_Z, int a_A, int a_M, bool bind ) {

    return( readTarget( a_lib_name, G4GIDI_Misc_Z_A_m_ToName( a_Z, a_A, a_M ), bind ) );
}

/* *********************************************************************************************************//**
 * Returns the protare specified by *a_targetName* and the projectile for *this*.
 *
 * @param a_lib_name            [in]    The evaluation. Call be an empty std::string.
 * @param a_targetName          [in]    The target PoPs id.
 * @param a_bind                [in]    If *true*, read target is added to member *m_protares*.
 *
 * @return                              Boolean value.
 ***********************************************************************************************************/

G4GIDI_target *G4GIDI::readTarget( std::string const &a_lib_name, std::string const &a_targetName, bool a_bind ) {

    for( auto iter_protare = m_protares.cbegin( ); iter_protare != m_protares.cend( ); ++iter_protare ) {
        if( *(*iter_protare)->getName( ) == a_targetName ) return( nullptr );
    }

    GIDI::Construction::Settings construction( GIDI::Construction::ParseMode::all, GIDI::Construction::PhotoMode::nuclearOnly );
    construction.setGRIN_continuumGammas( true );

    for( auto mapIter = m_maps.begin( ); mapIter!= m_maps.end( ); ++mapIter ) {
        GIDI::Protare *GIDI_protare = (*mapIter)->protare( construction, G4GIDI_pops, m_projectile, a_targetName, "", a_lib_name );
        if( GIDI_protare != nullptr ) {
            LUPI::StatusMessageReporting smr;
            GIDI::Styles::TemperatureInfos temperatures = GIDI_protare->temperatures( );
            std::string label( temperatures[0].griddedCrossSection( ) );
            MCGIDI::Transporting::MC MC( G4GIDI_pops, m_projectile, &GIDI_protare->styles( ), label, GIDI::Transporting::DelayedNeutrons::off, 20 ); // FIXME: 20
            MC.setThrowOnError( false );
            MC.setSampleNonTransportingParticles( true );
            MCGIDI::DomainHash domainHash( 4000, 1e-8, 10 );
            std::set<size_t> reactionsToExclude;

            GIDI::Transporting::Particles particles;
            GIDI::Transporting::Particle neutron( PoPI::IDs::neutron );
            particles.add( neutron );

            GIDI::Styles::TemperatureInfos temperatures1;
            temperatures1.push_back( temperatures[0] );
            MCGIDI::Protare *MCGIDI_protare = MCGIDI::protareFromGIDIProtare( smr, *GIDI_protare, G4GIDI_pops, MC, particles, domainHash, 
                    temperatures1, reactionsToExclude );
            //if( !smr.isOk( ) ) throw LUPI::Exception( smr.constructFullMessage( "G4GIDI::readTarget:" ) );

            G4GIDI_target *protare = new G4GIDI_target( G4GIDI_pops, domainHash, *GIDI_protare, MCGIDI_protare );
            delete GIDI_protare;
            if( a_bind ) m_protares.push_back( protare );
            return( protare );
        }
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Determines target name from *a_Z*, *a_A* and *a_M* and calls **getAlreadyReadTarget** with target's name, and returns
 * its return value.
 *
 * @param a_Z                   [in]    The atomic number of the target.
 * @param a_A                   [in]    The mass number of the taret.
 * @param a_M                   [in]    The meta-stable index of the target.
 *
 * @return                              Pointer to an **G4GIDI_target** instance of **nullptr**.
 ***********************************************************************************************************/

G4GIDI_target *G4GIDI::getAlreadyReadTarget( G4int a_Z, G4int a_A, G4int a_M ) {

    return( getAlreadyReadTarget( G4GIDI_Misc_Z_A_m_ToName( a_Z, a_A, a_M ) ) );
}

/* *********************************************************************************************************//**
 * Returns the pointer to the **G4GIDI_target** with name *a_targetName* in member *m_protares* or **nullptr** if one
 * does not exists in *m_protares*.
 *
 * @param a_targetName          [in]    The name of the target.
 *
 * @return                              Pointer to an **G4GIDI_target** instance of **nullptr**.
 ***********************************************************************************************************/

G4GIDI_target *G4GIDI::getAlreadyReadTarget( std::string const &a_targetName ) {

    for( auto iter_protare = m_protares.cbegin( ); iter_protare != m_protares.cend( ); ++iter_protare ) {
        if( *(*iter_protare)->getName( ) == a_targetName ) return( *iter_protare );
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * If *a_target* is in member *m_protares*, removed it from *m_protares*, delete it and return 0. Otherwise,
 * do nothing and return 1.
 *
 * @param a_target              [in]    The evaluation. Call be an empty std::string.
 *
 * @return                              Integer value.
 ***********************************************************************************************************/

G4int G4GIDI::freeTarget( G4GIDI_target *a_target ) {

    for( auto iter_protare = m_protares.cbegin( ); iter_protare != m_protares.cend( ); ++iter_protare ) {
        if( *iter_protare == a_target ) {
            m_protares.erase( iter_protare );
            delete a_target;
            return( 0 );
        }
    }
    return( 1 );
}

/* *********************************************************************************************************//**
 * Determines target name from *a_Z*, *a_A* and *a_M* and calls **freeTarget** with target's name, and returns
 * its return value.
 *
 * @param a_Z                   [in]    The atomic number of the target.
 * @param a_A                   [in]    The mass number of the taret.
 * @param a_M                   [in]    The meta-stable index of the target.
 *
 * @return                              Integer value.
 ***********************************************************************************************************/

G4int G4GIDI::freeTarget( G4int a_Z, G4int a_A, G4int a_M ) {

    return( freeTarget( G4GIDI_Misc_Z_A_m_ToName( a_Z, a_A, a_M ) ) );
}

/* *********************************************************************************************************//**
 * If target with name *targetSymbol* is in member *m_protares*, call **freeTarget** with target's pointer and return **freeTarget**'s return value.
 * Otherwise, return 1.
 *
 * @return                              Integer value.
 ***********************************************************************************************************/

G4int G4GIDI::freeTarget( std::string const &targetSymbol ) {

    for( auto iter_protare = m_protares.cbegin( ); iter_protare != m_protares.cend( ); ++iter_protare ) {
        if( *(*iter_protare)->getName( ) == targetSymbol ) return( freeTarget( *iter_protare ) );
    }
    return( 1 );
}

/* *********************************************************************************************************//**
 * If target with name *targetSymbol* is in member *m_protares*, call **freeTarget** with target's pointer and return **freeTarget**'s return value.
 * Otherwise, return 1.
 *
 * @return                              Integer value.
 ***********************************************************************************************************/

std::vector<std::string> *G4GIDI::getListOfReadTargetsNames( void ) {

    std::vector<std::string> *listOfTargets = new std::vector<std::string>( );

    for( auto iter_protare = m_protares.cbegin( ); iter_protare != m_protares.cend( ); ++iter_protare ) {
        listOfTargets->push_back( *(*iter_protare)->getName( ) );
    }

    return( listOfTargets );
}
