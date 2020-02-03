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
//
// 
// Laurent Garnier  27th October 2011

#if defined (G4VIS_BUILD_OPENGLQT_DRIVER) || defined (G4VIS_USE_OPENGLQT)

#ifndef G4OPENGLSTOREDQTSCENEHANDLER_HH
#define G4OPENGLSTOREDQTSCENEHANDLER_HH

#include "G4OpenGLStoredSceneHandler.hh"

class G4OpenGLStoredQtSceneHandler: public G4OpenGLStoredSceneHandler {

public:

  G4OpenGLStoredQtSceneHandler (G4VGraphicsSystem& system, const G4String& name = "");
  virtual ~G4OpenGLStoredQtSceneHandler ();

  // Two virtual functions for extra processing in a sub-class, for
  // example, to make a display tree.  They are to return true if the
  // visible object uses gl commands for drawing.  This is
  // predominantly true; a notable exception is Qt text.  In that
  // case, a display list does not need to be created; all relevant
  // information is assumed to be stored in the PO/TOList.
  G4bool ExtraPOProcessing(const G4Visible&, size_t currentPOListIndex);
  G4bool ExtraTOProcessing(const G4Visible&, size_t currentTOListIndex);

  void ClearStore ();
  void ClearTransientStore ();
  void SetScene(G4Scene*);
};

#endif

#endif
