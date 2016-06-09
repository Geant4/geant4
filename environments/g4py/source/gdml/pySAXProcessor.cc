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
// $Id: pySAXProcessor.cc,v 1.1 2007/11/13 10:03:23 kmura Exp $
// $Name: geant4-09-01 $
// ====================================================================
//   pySAXProcessor.cc
//
//   based on 2.10.0
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "Saxana/SAXProcessor.h"
#include "Saxana/ProcessingConfigurator.h"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_SAXProcessor()
{
  class_<SAXProcessor>("SAXProcessor", "SAX processor class")
    // ---
    .def("Initialize",   &SAXProcessor::Initialize)
    .def("Configure",    &SAXProcessor::Configure)
    .def("Run",          &SAXProcessor::Run)
    .def("Finalize",     &SAXProcessor::Finalize)
    ;
}

