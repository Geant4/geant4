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
#ifndef HEPVis_SoGL2PSAction_h
#define HEPVis_SoGL2PSAction_h 

#include <Inventor/actions/SoGLRenderAction.h>

/**
 *  SoGL2PSAction inherits Inventor/SoGLRenderAction.
 * It uses the gl2ps software to produce vector PostScript of a scene graph.
 *  See applications/DetectorTreeKit, Polyhedron for examples.
 */

#define SoGL2PSAction Geant4_SoGL2PSAction

class SoGL2PSAction : public SoGLRenderAction {
  SO_ACTION_HEADER(SoGL2PSAction);
public:
  SoGL2PSAction(const SbViewportRegion&);
  void setFileName(const char*);
  void enableFileWriting();
  void disableFileWriting();
  SbBool fileWritingEnabled() const;
  SbBool addBitmap(int,int,float=0,float=0,float=0,float=0);
  void beginViewport();
  void endViewport();
public: /*SoINTERNAL*/
  static void initClass();
protected:
  virtual void beginTraversal(SoNode*);
private:
  void gl2psBegin();
private:
  SbString fFileName;
  FILE* fFile;
};

#endif

