// $Id: G4PyCoutDestination.hh,v 1.1 2006-04-25 08:13:51 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   G4PyCoutDestination.hh
//
//                                         2006 Q
// ====================================================================
#ifndef G4PY_COUT_DESTINATION_H
#define G4PY_COUT_DESTINATION_H

#include "G4coutDestination.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class G4PyCoutDestination : public G4coutDestination {

public:
  G4PyCoutDestination();
  ~G4PyCoutDestination();

  virtual G4int ReceiveG4cout(G4String coutString);
  virtual G4int ReceiveG4cerr(G4String cerrString);

};

#endif
