// $Id: pyG4UImanager.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4UImanager.cc
//
//   G4UImanager class is pure singleton, so it cannnot be exposed
//   from BPL. Functionality of G4UImanager is exposed in global
//   name space via wrappers.
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UImanager.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
namespace pyG4UImanager {

//////////////////////////////////////////////
G4int ApplyUICommand_1(const G4String& cmdstr)
//////////////////////////////////////////////
{
  G4UImanager* UImgr= G4UImanager::GetUIpointer();
  G4int returnVal= UImgr-> ApplyCommand(cmdstr);
  if( returnVal == fCommandSucceeded ) return returnVal;

  G4int paramIndex= returnVal % 100;
  G4int commandStatus= returnVal - paramIndex;
 
  switch(commandStatus) {
    case fCommandSucceeded:
      break;

    case fCommandNotFound:
      G4cout << "command <" << UImgr-> SolveAlias(cmdstr)
	     << "> not found" << G4endl;
      break;

    case fIllegalApplicationState:
      G4cout << "illegal application state -- command refused" 
	     << G4endl;
      break;

    case fParameterOutOfRange:
      break;

    case fParameterOutOfCandidates:
      G4cout << "Parameter is out of candidate list (index "
	     << paramIndex << ")"
	     << G4endl;
      break;

    case fParameterUnreadable:
      G4cout << "Parameter is wrong type and/or is not omittable (index " 
	     << paramIndex << ")" << G4endl;
      break;

    case fAliasNotFound:
      break;

    default:
      G4cout << "command refused (" << commandStatus << ")" << G4endl;
      break;
  }

  return returnVal;
}

/////////////////////////////////////////////////
G4int ApplyUICommand_2(const std::string& cmdstr)
/////////////////////////////////////////////////
{
  return ApplyUICommand_1(cmdstr);
}


////////////////////////////////////////////////
std::string GetCurrentValues(const char* cmdstr)
////////////////////////////////////////////////
{
  G4UImanager* UImgr= G4UImanager::GetUIpointer();
  return UImgr-> GetCurrentValues(cmdstr);  
}

};

using namespace pyG4UImanager;

// ====================================================================
// module definition
// ====================================================================
void export_G4UImanager()
{
  def("ApplyUICommand",    ApplyUICommand_1);
  def("ApplyUICommand",    ApplyUICommand_2);
  def("GetCurrentValues",  GetCurrentValues);
}

