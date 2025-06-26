/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "HAPI.hpp"

namespace HAPI {

/*
=========================================================
 *
 * @return
 */
Text::Text() :
        m_text( "" ) {

}
/*
=========================================================
 *
 * @param a_text text string
 * @return
 */
Text::Text( std::string const a_text ) :
        m_text( a_text ) {

}
/*
=========================================================
*/
Text::~Text( ) {

}

}

