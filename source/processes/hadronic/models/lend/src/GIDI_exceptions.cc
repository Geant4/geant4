/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <sstream>
#include <algorithm>

#include "GIDI.hpp"

namespace GIDI {

/*! \class Exception
 * Exception class for all GIDI exceptions thrown by GIDI functions.
 */

/* *********************************************************************************************************//**
 * @param a_message         [in]     The message that the function what() will return.
 ***********************************************************************************************************/

Exception::Exception( std::string const & a_message ) :
        std::runtime_error( a_message ) {

}

}               // End namespace GIDI.
