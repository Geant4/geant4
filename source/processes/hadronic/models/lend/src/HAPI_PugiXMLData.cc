/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>

#include "HAPI.hpp"

#ifdef HAPI_USE_PUGIXML
namespace HAPI {

/*
=========================================================
 *
 * @return
 */
PugiXMLData::PugiXMLData() :
        m_node( pugi::xml_node() ),
        m_length( -1 ) {

}
/*
=========================================================
 *
 * @param a_node
 * @return
 */
PugiXMLData::PugiXMLData( pugi::xml_node a_node ) :
        m_node( a_node ),
        m_length( -1 ) {

}
/*
=========================================================
*/
PugiXMLData::~PugiXMLData( ) {

}

int PugiXMLData::length( ) const {

    if (m_length == -1)
        throw "Can't access length until data is read!";
    return m_length;

}

void PugiXMLData::getDoubles(nf_Buffer<double> &buffer)
{
    int64_t numberConverted;
    char *endCharacter;
    char const *text = m_node.text( ).get( );
    double *dValues = nfu_stringToListOfDoubles( NULL, text, ' ', &numberConverted, &endCharacter, 0 );
    if (dValues == NULL) throw "dValues = NULL";
    if (*endCharacter != 0) throw "bad values string";
    m_length = numberConverted;

    buffer.resize(m_length);
    for (int i=0; i<m_length; i++)
        buffer[i] = dValues[i];
    free( dValues );
}

void PugiXMLData::getInts(nf_Buffer<int> &buffer)
{
    int64_t numberConverted;
    char *endCharacter;
    char const *text = m_node.text( ).get( );
    int *iValues = nfu_stringToListOfInt32s( NULL, text, ' ', &numberConverted, &endCharacter );
    if (iValues == NULL) throw "dValues = NULL";
    if (*endCharacter != 0) throw "bad values string";
    m_length = numberConverted;

    buffer.resize(m_length);
    for (int i=0; i<m_length; i++)
        buffer[i] = iValues[i];
    free( iValues );
}

}
#endif
