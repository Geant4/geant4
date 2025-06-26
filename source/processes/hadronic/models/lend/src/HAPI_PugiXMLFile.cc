/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "HAPI.hpp"

#ifdef HAPI_USE_PUGIXML
namespace HAPI {

/*
=========================================================
 *
 * @return
 */
PugiXMLFile::PugiXMLFile() :
        m_name( "" ) {

}
/*
=========================================================
 *
 * @param filename
 * @return
 */
PugiXMLFile::PugiXMLFile(char const *filename, std::string const &a_callingFunctionName) :
        m_name( filename ) {

    pugi::xml_parse_result result = m_doc.load_file(filename);
    if (result.status != pugi::status_ok) {
        throw std::runtime_error( "ERROR from PugiXMLFile::PugiXMLFile via " + a_callingFunctionName + " for file '" + filename + "': " + result.description() );
    }

}
/*
=========================================================
*/
PugiXMLFile::~PugiXMLFile( ) {

}
/*
============================================================
===================== get child element ====================
============================================================
 *
 * @return
 */
Node PugiXMLFile::child(char const *a_name) {

    // Only one child element allowed in XML file
    if (a_name == m_doc.first_child( ).name())
        return Node(new PugiXMLNode( m_doc.first_child( ) ));
    else
        return Node(new PugiXMLNode( ));
}
/*
=========================================================
*/
Node PugiXMLFile::first_child() {

    return Node(new PugiXMLNode( m_doc.first_child( ) ));

}
/*
=========================================================
*/
std::string PugiXMLFile::name() const {

    return m_name;

}

}
#endif
