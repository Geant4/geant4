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
//
// $Id: G4Wt.hh,v 1.6 2010-05-20 07:01:03 lgarnier Exp $
//
//  To unify Wt event treatment between 
// G4/interfaces Wt sessions and G4/visualizations Wt drivers.
// L. Garnier
#ifndef G4WT_HH
#define G4WT_HH

#if defined(G4INTY_BUILD_WT) || defined(G4INTY_USE_WT)

#include "G4VInteractorManager.hh"

// Class description :
//
//  G4Wt : a singleton to handle GUI sessions and visualization 
// drivers built over Wt. It permits to have one Wt main loop for 
// the whole application. The Wt toolkit is inited in the 
// constructor. It is done once for the whole application.
//
// Class description - end :

#include <Wt/WServer>

class G4Wt : public G4VInteractorManager {
public:
  static G4Wt* getInstance();
  static G4Wt* getInstance(int,char**,char*);
  Wt::WServer getServer();
  G4bool Inited();
  void* GetEvent();
  void FlushAndWaitExecution();
  virtual ~G4Wt();                     
private:
  G4Wt (int,char**,char*);                     
  static G4Wt* instance; // Pointer to single instance.
  Wt::WServer *wServer;
  int    argn;
  char** args;
};

#endif //HAS_WT

#endif
