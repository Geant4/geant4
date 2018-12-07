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
// F.W. Jones 05012018

#ifndef G4SOQT_HH
#define G4SOQT_HH

#if defined(G4VIS_BUILD_OIQT_DRIVER)
//#if defined(G4INTY_BUILD_QT) || defined(G4INTY_USE_QT)

#include "G4VInteractorManager.hh"

class QWidget;
class QString;

// Class description :
//
// G4SoQt : a singleton to handle Open Inventor visualization 
// drivers built over Qt. It permits to have one Qt main loop for 
// the whole application. The SoQt toolkit is inited in the 
// constructor. It is done once for the whole application.
//
// Class description - end :

class G4SoQt : public G4VInteractorManager {
public:
  static G4SoQt* getInstance();
  // FWJ no command line args for time being
  //  static G4SoQt* getInstance(int, char**, char*);
  G4bool Inited();
  void* GetEvent();
  void FlushAndWaitExecution();
  void SecondaryLoop();
  virtual ~G4SoQt();                     
  bool IsExternalApp();

private:
  G4SoQt(const G4SoQt&);
  G4SoQt();
  //  G4SoQt(int, char**, char*);                     
  G4SoQt& operator=(const G4SoQt&);
  static G4SoQt* instance; // Pointer to single instance.
  //  int    argn;
  //  char** args;
  bool externalApp;
};

#endif

#endif
