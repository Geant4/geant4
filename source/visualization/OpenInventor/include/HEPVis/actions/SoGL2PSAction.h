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

#include <tools/gl2ps_def.h>

#include <string>

/**
 *  SoGL2PSAction inherits Inventor/SoGLRenderAction.
 * It uses the gl2ps software to produce vector PostScript of a scene graph.
 */

#define SoGL2PSAction Geant4_SoGL2PSAction

class SoGL2PSAction : public SoGLRenderAction {
  SO_ACTION_HEADER(SoGL2PSAction);
public:
  SoGL2PSAction(const SbViewportRegion&);
  virtual ~SoGL2PSAction();
public: /*SoINTERNAL*/
  static void initClass();
public:  
  void setFileName(const std::string&);
  void setTitleAndProducer(const std::string&,const std::string&);
  void setExportImageFormat_PS();
  void setExportImageFormat_EPS();
  void setExportImageFormat_TEX();
  void setExportImageFormat_PDF();
  void setExportImageFormat_SVG();
  void setExportImageFormat_PGF();
  bool enableFileWriting();
  void disableFileWriting();
  bool addBitmap(int,int,float=0,float=0,float=0,float=0);
protected:
  virtual void beginTraversal(SoNode*);
private:  
  bool openFile();
  void closeFile();
  bool beginPage(int,int,int,int);
  bool endPage();
  tools_GL2PScontextPointer fContext;
  FILE* fFile;
  std::string fFileName;
  std::string fTitle;
  std::string fProducer;
  unsigned int fFormat;
};

#endif

