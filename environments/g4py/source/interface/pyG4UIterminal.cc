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
// $Id: pyG4UIterminal.cc,v 1.8 2008-03-13 07:32:18 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UIterminal.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

using namespace boost::python;

static G4UIterminal* session= 0;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4UIterminal {

/////////////////////
void StartUISession()
/////////////////////
{
  if (session == 0 ) {
    G4UItcsh* tcsh= new 
      G4UItcsh("[40;01;33mg4py[40;31m(%s)[40;36m[%/][00;30m:");

#if G4VERSION_NUMBER >= 900
    session= new G4UIterminal(tcsh, false);
#else
    session= new G4UIterminal(tcsh);
#endif
    tcsh-> SetLsColor(BLUE,RED);
  }

  session-> SessionStart();
}

};

using namespace pyG4UIterminal;

// ====================================================================
// module definition
// ====================================================================
void export_G4UIterminal()
{
  def("StartUISession", StartUISession);
}

