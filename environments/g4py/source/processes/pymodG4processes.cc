// $Id: pymodG4processes.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4processes.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4ProcessManager();
void export_G4ProcessTable();
void export_G4VProcess();
void export_G4ProcVector();
void export_G4ProcessType();

BOOST_PYTHON_MODULE(G4processes)
{
  export_G4ProcessManager();
  export_G4ProcessTable();
  export_G4VProcess();
  export_G4ProcVector();
  export_G4ProcessType();
}

