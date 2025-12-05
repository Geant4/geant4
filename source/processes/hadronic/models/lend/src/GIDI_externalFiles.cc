/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "GIDI.hpp"

namespace GIDI {

namespace ExternalFiles {

/*!
 *  loop over external files, if any represent binary store then open up HAPI::DataManager for that file
 */
void Suite::registerBinaryFiles( LUPI_maybeUnused std::string const &a_parentDir, LUPI_maybeUnused SetupInfo &a_setupInfo ) {

#ifdef HAPI_USE_HDF5
    if (this->has( "HDF" )) {
        a_setupInfo.m_protare->setDataManager( new HAPI::HDFDataManager( a_parentDir + "/" + this->get<ExternalFile>( "HDF" )->path( ) ) );
    }
#endif
}

}

}
