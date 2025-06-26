/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "MCGIDI.hpp"

namespace MCGIDI {

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void URR_protareInfo::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_CAST( m_inURR, a_buffer, a_mode, bool );
    DATA_MEMBER_DOUBLE( m_rng_Value, a_buffer, a_mode );
}

/* *********************************************************************************************************//**
 * URR_protareInfos constructor.
 *
 * @param a_protares            [in]    The list of protares to be check for URR data. Each protare with URR data add to *a_URR_protareInfos*.
 ***********************************************************************************************************/

LUPI_HOST URR_protareInfos::URR_protareInfos( Vector<Protare *> &a_protares ) {

    setup( a_protares );
}

/* *********************************************************************************************************//**
 * URR_protareInfos setup.
 *
 * @param a_protares            [in]    The list of protares to be check for URR data. Each protare with URR data add to *a_URR_protareInfos*.
 ***********************************************************************************************************/

LUPI_HOST void URR_protareInfos::setup( Vector<Protare *> &a_protares ) {

    std::vector<URR_protareInfo> URR_protareInfo_1;

    for( std::size_t i1 = 0; i1 < a_protares.size( ); ++i1 ) {
        Protare *protare = a_protares[i1];

        for( std::size_t i2 = 0; i2 < protare->numberOfProtares( ); ++i2 ) {
            ProtareSingle *protareSingle = const_cast<ProtareSingle *>( protare->protare( i2 ) );

            if( protareSingle->hasURR_probabilityTables( ) ) {
                protareSingle->URR_index( URR_protareInfo_1.size( ) );
                URR_protareInfo_1.push_back( URR_protareInfo( ) );
            }
        }
    }

    m_URR_protareInfos.reserve( URR_protareInfo_1.size( ) );
    m_URR_protareInfos.clear( );
    for( std::size_t i1 = 0; i1 < URR_protareInfo_1.size( ); ++i1 ) m_URR_protareInfos.push_back( URR_protareInfo_1[i1] );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void URR_protareInfos::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    std::size_t vectorSize = m_URR_protareInfos.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_URR_protareInfos.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_URR_protareInfos.internalSize();

    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_URR_protareInfos[vectorIndex].serialize( a_buffer, a_mode );
    }
}

/*! \class ACE_URR_probabilityTable
 * Class to store ACE URR probability table at one projectile energy for one type of reaction (e.g., total, elastic).
 */

/* *********************************************************************************************************//**
 * Simple constructor needed for broadcasting.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ACE_URR_probabilityTable::ACE_URR_probabilityTable( ) :
        m_energy( 0.0 ) {

}

/* *********************************************************************************************************//**
 * @param a_energy              [in]    The projectile energy where the data are specified.
 * @param a_propabilities       [in]    The probability for each cross section.
 * @param a_crossSection        [in]    The cross section for each probability.
 ***********************************************************************************************************/

LUPI_HOST ACE_URR_probabilityTable::ACE_URR_probabilityTable( double a_energy, std::vector<double> const &a_propabilities, 
                std::vector<double> const &a_crossSection ) :
        m_energy( a_energy ),
        m_propabilities( a_propabilities ),
        m_crossSections( a_crossSection ) {

    double sum = 0.0;
    for( std::size_t index = 0; index < m_propabilities.size( ); ++index ) {
        sum += m_propabilities[index];
        m_propabilities[index] = sum;
    }
    m_propabilities[m_propabilities.size( )-1] = 1.0;
}

/* *********************************************************************************************************//**
 * Simple constructor needed for broadcasting.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ACE_URR_probabilityTable::~ACE_URR_probabilityTable( ) {

}

/* *********************************************************************************************************//**
 * Returns the cross section corresponding to the probability *a_rng_Value*.
 *
 * @param a_rng_Value           [in]    A random number in the range [0,1).
 *
 * @return                              The cross section associated with probability a_rng_Value;
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ACE_URR_probabilityTable::sample( double a_rng_Value ) {

    int index = binarySearchVector( a_rng_Value, m_propabilities, true );
    if( m_propabilities[index] < a_rng_Value ) ++index;
    return( m_crossSections[index] );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ACE_URR_probabilityTable::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_DOUBLE( m_energy, a_buffer, a_mode  );
    DATA_MEMBER_VECTOR_DOUBLE( m_propabilities, a_buffer, a_mode  );
    DATA_MEMBER_VECTOR_DOUBLE( m_crossSections, a_buffer, a_mode  );
}

/*! \class ACE_URR_probabilityTables
 * Class to store ACE URR probability tables at a list of projectile energies for one type of reaction (e.g., total, elastic).
 */

/* *********************************************************************************************************//**
 * Simple constructor needed for broadcasting.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ACE_URR_probabilityTables::ACE_URR_probabilityTables( ) {

}

/* *********************************************************************************************************//**
 * @param a_capacity            [in]    The number of energy slots to reverse.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ACE_URR_probabilityTables::ACE_URR_probabilityTables( std::size_t a_capacity ) {

    m_energies.reserve( a_capacity );
    m_ACE_URR_probabilityTables.reserve( a_capacity );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ACE_URR_probabilityTables::~ACE_URR_probabilityTables( ) {

    for( auto iter = m_ACE_URR_probabilityTables.begin( ); iter != m_ACE_URR_probabilityTables.end( ); ++iter ) delete (*iter);
}

/* *********************************************************************************************************//**
 * Calls reserve for m_energies and m_ACE_URR_probabilityTables with the value *a_capacity*.
 *
 * @param a_capacity                    [in]    The size of the space to reserve.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ACE_URR_probabilityTables::reserve( std::size_t a_capacity ) {

    m_energies.reserve( a_capacity );
    m_ACE_URR_probabilityTables.reserve( a_capacity );
}

/* *********************************************************************************************************//**
 * Adds *a_ACE_URR_probabilityTable* to the end of *this*.
 *
 * @param a_ACE_URR_probabilityTable    [in]    **ACE_URR_probabilityTable** instance to add.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ACE_URR_probabilityTables::push_back( ACE_URR_probabilityTable *a_ACE_URR_probabilityTable ) {

    if( m_energies.size( ) == capacity( ) ) LUPI_THROW( "ACE_URR_probabilityTables::addEnergyData: adding too many ACE_URR_probabilityTables." );
    m_energies.push_back( a_ACE_URR_probabilityTable->energy( ) );
    m_ACE_URR_probabilityTables.push_back( a_ACE_URR_probabilityTable );
}

/* *********************************************************************************************************//**
 * Returns the cross section corresponding to the probability *a_rng_Value*.
 *
 * @param a_energy              [in]    The incident projectiles energy.
 * @param a_rng_Value           [in]    A random number in the range [0,1).
 *
 * @return                              The cross section associated with probability a_rng_Value;
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ACE_URR_probabilityTables::sample( double a_energy, double a_rng_Value ) {

    int index = binarySearchVector( a_energy, m_energies, true );

    std::size_t index_t = (std::size_t) index;
    if( index_t < m_energies.size( ) - 1 ) {
        if( 0.5 * ( m_energies[index_t] + m_energies[index_t+1] ) < a_energy ) ++index_t;     // Find closest energy.
    }

    return( m_ACE_URR_probabilityTables[index_t]->sample( a_rng_Value ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ACE_URR_probabilityTables::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_VECTOR_DOUBLE( m_energies, a_buffer, a_mode );

    std::size_t vectorSize = m_energies.size( );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_ACE_URR_probabilityTables.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_ACE_URR_probabilityTables.internalSize( );

    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        ACE_URR_probabilityTable *ACE_URR_probabilityTable1 = m_ACE_URR_probabilityTables[vectorIndex];
        if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
            if( a_buffer.m_placement != nullptr ) {
                ACE_URR_probabilityTable1 = new(a_buffer.m_placement) ACE_URR_probabilityTable;
                a_buffer.incrementPlacement( sizeof( ACE_URR_probabilityTable ) ); }
            else {
                ACE_URR_probabilityTable1 = new ACE_URR_probabilityTable;
            }
            m_ACE_URR_probabilityTables[vectorIndex] = ACE_URR_probabilityTable1;
        }
        if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof( ACE_URR_probabilityTable ) );
        }
        ACE_URR_probabilityTable1->serialize( a_buffer, a_mode );
    }
}

/* *********************************************************************************************************//**
 * This method serializes a **ACE_URR_probabilityTables** instance pointed to by *a_ACE_URR_probabilityTables*.
 *
 * @param a_protare             [in]    The GIDI::Protare whose data is to be used to construct *this*.
 * @param a_settings            [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_setupInfo           [in]    Used internally when constructing a Protare to pass information to other constructors.
 *
 * @return                              A pointer to the converted **ACE_URR_probabilityTables** instance.
 ***********************************************************************************************************/

LUPI_HOST void convertACE_URR_probabilityTablesFromGIDI( GIDI::ProtareSingle const &a_protare, Transporting::MC &a_settings, SetupInfo &a_setupInfo ) {

    if( ( a_settings.crossSectionLookupMode( ) == Transporting::LookupMode::Data1d::continuousEnergy ) 
            && ( a_settings._URR_mode( ) == Transporting::URR_mode::ACE_URR_probabilityTables ) ) {

        int64_t numberConverted;
        char *endCharacter;

        for( auto iter = a_protare.ACE_URR_probabilityTables( ).begin( ); iter != a_protare.ACE_URR_probabilityTables( ).end( ); ++iter ) {
            bool needToInitialize( true );
            std::map<int, std::string> columnNames;
            ACE_URR_probabilityTablesFromGIDI *ACE_URR_probabilityTablesFromGIDI1 = new ACE_URR_probabilityTablesFromGIDI( );
            GIDI::ACE_URR::ProbabilityTable *form = dynamic_cast<GIDI::ACE_URR::ProbabilityTable *>( *iter );
            GIDI::ACE_URR::ProbabilityTable::Forms &incidentEnergies = form->forms( );

            for( auto incidentEnergyIter = incidentEnergies.begin( ); incidentEnergyIter != incidentEnergies.end( ); ++incidentEnergyIter ) {
                GIDI::ACE_URR::IncidentEnergy *incidentEnergy = *incidentEnergyIter;
                GIDI::Table::Table const &table = incidentEnergy->table( );
                int numberOfRows = table.rows( );
                int numberOfColumns = table.columns( );

                int columnIndex = 0;
                for( auto columnHeaderIter = table.columnHeaders( ).begin( ); columnHeaderIter != table.columnHeaders( ).end( ); ++columnHeaderIter ) {
                    GIDI::Table::Column const *columnHeader = dynamic_cast<GIDI::Table::Column *>( *columnHeaderIter );
                    if( needToInitialize && columnIndex > 0 ) {
                        ACE_URR_probabilityTablesFromGIDI1->m_ACE_URR_probabilityTables[columnHeader->name( )] = 
                                new ACE_URR_probabilityTables( incidentEnergies.size( ) );
                        columnNames[columnIndex] = columnHeader->name( );
                    }
                    ++columnIndex;
                }
                needToInitialize = false;

                GIDI::Table::Data const &data = table.data( );

                std::string const &body = data.body( );
                char const *text = body.c_str( );
                double *dValues = nfu_stringToListOfDoubles( NULL, text, data.sep( ).c_str( )[0], &numberConverted, &endCharacter, 0 );
                if( dValues == nullptr ) throw GIDI::Exception( "convertACE_URR_probabilityTablesFromGIDI: nfu_stringToListOfDoubles failed." );

                std::vector<std::vector<double> > columns( numberOfColumns );
                for( columnIndex = 0; columnIndex < numberOfColumns; ++columnIndex ) {
                    columns[columnIndex].reserve( numberOfRows );
                    for( int rowIndex = 0; rowIndex < numberOfRows; ++rowIndex ) columns[columnIndex].push_back( dValues[rowIndex*numberOfColumns+columnIndex] );
                }
                free( dValues );

                for( columnIndex = 1; columnIndex < numberOfColumns; ++columnIndex ) {
                    ACE_URR_probabilityTablesFromGIDI1->m_ACE_URR_probabilityTables[columnNames[columnIndex]]->push_back( 
                            new ACE_URR_probabilityTable( incidentEnergy->value( ), columns[0], columns[columnIndex] ) );
                }
            }
            a_setupInfo.m_ACE_URR_probabilityTablesFromGIDI[form->label()] = ACE_URR_probabilityTablesFromGIDI1;
        }
    }
}

/* *********************************************************************************************************//**
 * This method serializes a **Transporting::URR_mode** value.
 *
 * @param a_URR_mode            [in]    The inputted Transporting::URR_mode value.
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 *
 * @return                              The Transporting::URR_mode value.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Transporting::URR_mode serializeURR_mode( Transporting::URR_mode a_URR_mode, LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    int type = 0;
    switch( a_URR_mode ) {
    case Transporting::URR_mode::none :
        break;
    case Transporting::URR_mode::pdfs :
        type = 1;
        break;
    case Transporting::URR_mode::ACE_URR_probabilityTables :
        type = 2;
        break;
    }
    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( type == 0 ) return( Transporting::URR_mode::none );
    if( type == 1 ) return( Transporting::URR_mode::pdfs );
    return( Transporting::URR_mode::ACE_URR_probabilityTables );
}

/* *********************************************************************************************************//**
 * This method serializes a **ACE_URR_probabilityTables** instance pointed to by *a_ACE_URR_probabilityTables*.
 *
 * @param a_ACE_URR_probabilityTables   [in]    Specifies the action of this method.
 * @param a_buffer                      [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                        [in]    Specifies the action of this method.
 *
 * @return                              A pointer to the serialized **ACE_URR_probabilityTables** instance.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ACE_URR_probabilityTables *serializeACE_URR_probabilityTables( ACE_URR_probabilityTables *a_ACE_URR_probabilityTables, 
                LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    int type = 0;
    if( a_ACE_URR_probabilityTables != nullptr ) type = 1;
    DATA_MEMBER_INT( type, a_buffer, a_mode );
    if( type == 0 ) return( nullptr );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        if( a_buffer.m_placement != nullptr ) {
            a_ACE_URR_probabilityTables = new(a_buffer.m_placement) ACE_URR_probabilityTables;
            a_buffer.incrementPlacement( sizeof( ACE_URR_probabilityTables ) ); }
        else {
            a_ACE_URR_probabilityTables = new ACE_URR_probabilityTables;
        }
    }

    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.incrementPlacement( sizeof( ACE_URR_probabilityTables ) );

    a_ACE_URR_probabilityTables->serialize( a_buffer, a_mode );

    return( a_ACE_URR_probabilityTables );
}

/*! \class ACE_URR_probabilityTablesFromGIDI
 * Class to store temporary ACE URR probability table data.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST ACE_URR_probabilityTablesFromGIDI::ACE_URR_probabilityTablesFromGIDI( ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST ACE_URR_probabilityTablesFromGIDI::~ACE_URR_probabilityTablesFromGIDI( ) {

    for( auto iter = m_ACE_URR_probabilityTables.begin( ); iter != m_ACE_URR_probabilityTables.end( ); ++iter ) delete (*iter).second;

}

}       // End namespace MCGIDI.
