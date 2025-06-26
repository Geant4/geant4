/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <string.h>
#include <limits.h> 

#include "GIDI.hpp"

namespace GIDI {

namespace Construction {

/*! \class Settings
 * This class is used to pass user parameters to various constructors.
 *
 * The main use is to limit the type of data read in via the **a_parseMode** argument (see enum ParseMode).
*/

/* *********************************************************************************************************//**
 * @param a_parseMode           [in]    Instructs the parses on which data to parse.
 * @param a_photoMode           [in]    Instructs the parses if photo atomic and/or photoicnuclear protares are to be included.
 ***********************************************************************************************************/

Settings::Settings( ParseMode a_parseMode, PhotoMode a_photoMode ) :
        m_parseMode( a_parseMode ),
        m_photoMode( a_photoMode ),
        m_useSystem_strtod( 0 ),
        m_lazyParsing( true ),
        m_decayPositronium( true ),
        m_usePhotoAtomicIncoherentDoppler( false ),
        m_fissionResiduals( FissionResiduals::none ),
        m_GRIN_continuumGammas( false ) {

}

/* *********************************************************************************************************//**
 * Copy constructor.
 *
 * @param a_settings            [in]    The **Settings** instance to copy.
 ***********************************************************************************************************/

Settings::Settings( Settings const &a_settings ) :
        m_parseMode( a_settings.parseMode( ) ),
        m_photoMode( a_settings.photoMode( ) ),
        m_useSystem_strtod( a_settings.useSystem_strtod( ) ),
        m_lazyParsing( a_settings.lazyParsing( ) ),
        m_decayPositronium( a_settings.decayPositronium( ) ),
        m_usePhotoAtomicIncoherentDoppler( a_settings.usePhotoAtomicIncoherentDoppler( ) ),
        m_fissionResiduals( a_settings.fissionResiduals( ) ),
        m_GRIN_continuumGammas( false ) {

}

}

}
