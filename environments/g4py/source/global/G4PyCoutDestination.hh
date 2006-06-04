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
// $Id: G4PyCoutDestination.hh,v 1.2 2006-06-04 21:34:29 kmura Exp $
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
