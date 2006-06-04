//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: pyG4UIterminal.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UIterminal.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
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
      G4UItcsh("[40;01;33meccsim[40;31m(%s)[40;36m[%/][00;01;30m:");

    session= new G4UIterminal(tcsh);
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

