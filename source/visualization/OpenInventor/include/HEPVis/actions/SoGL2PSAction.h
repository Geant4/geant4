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

