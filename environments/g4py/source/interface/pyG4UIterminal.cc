// $Id: pyG4UIterminal.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
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

