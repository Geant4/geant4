//
// File name:     RadmonApplicationRunNumbering.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationRunNumbering.cc,v 1.1 2005-11-24 02:34:21 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplicationRunNumbering.hh"
#include "G4Run.hh"

void                                            RadmonApplicationRunNumbering :: OnBeginOfRun(const G4Run * run)
{
 if (!enable)
  return;
  
 G4cout << "RadmonApplicationRunNumbering::OnBeginOfRun: Begin of run " << run->GetRunID() << '.' << G4endl;
}
