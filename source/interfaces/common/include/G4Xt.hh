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
// $Id: G4Xt.hh 70601 2013-06-03 11:20:53Z gcosmo $
//
//  To unify X11 event treatment between 
// G4/interfaces Xt sessions and G4/visualizations Xt drivers.
// G.Barrand

#ifndef G4XT_HH
#define G4XT_HH

#if defined(G4INTY_BUILD_XT) || defined(G4INTY_USE_XT)

#include <X11/Intrinsic.h>

#include "G4VInteractorManager.hh"

// Class description :
//
//  G4Xt : a singleton to handle GUI sessions and visualization 
// drivers built over Xt. It permits to have one Xt main loop for 
// the whole application. The Xt toolkit is inited in the 
// constructor. It is done once for the whole application.
//
// Class description - end :

class G4Xt : public G4VInteractorManager {
public:
  static G4Xt* getInstance();
  static G4Xt* getInstance(int,char**,char*);
  void PutStringInResourceDatabase(char*);
  G4bool Inited();
  void* GetEvent();
  void FlushAndWaitExecution();
  virtual ~G4Xt();                     
private:
  G4Xt (const G4Xt&);
  G4Xt (int,char**,char*);
  G4Xt& operator= (const G4Xt&);
  static G4Xt* instance; // Pointer to single instance.
};

#endif //HAS_XT

#endif
