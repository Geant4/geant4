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
// $Id$
//
/// \file Analysis.hh
/// \brief Selection of the analysis technology

#ifndef Analysis_h
#define Analysis_h 1

// Default output format

#if ( ( !defined(TEST_ANALYSIS_CSV) ) && ( !defined(TEST_ANALYSIS_HBOOK) ) &&\
      ( !defined(TEST_ANALYSIS_ROOT) ) && ( !defined(TEST_ANALYSIS_XML) ) )
#define TEST_ANALYSIS_ROOT 1
#endif      

#ifdef TEST_ANALYSIS_ROOT
#include "g4root.hh"
#endif

#ifdef TEST_ANALYSIS_XML
#include "g4xml.hh"
#endif

#ifdef TEST_ANALYSIS_CSV
#include "g4csv.hh"
#endif

#ifdef TEST_ANALYSIS_HBOOK
#include "g4hbook.hh"
#endif

#endif
