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
// $Id: G4OpenGL2PSAction.cc,v 1.6 2010-11-03 16:40:34 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifdef G4VIS_BUILD_OPENGL_DRIVER
 #define G4VIS_BUILD_OPENGL_GL2PS 
#endif
#ifdef G4VIS_BUILD_OI_DRIVER
 #define G4VIS_BUILD_OPENGL_GL2PS 
#endif

#ifdef G4VIS_BUILD_OPENGL_GL2PS

#include "G4OpenGL2PSAction.hh"


G4OpenGL2PSAction::G4OpenGL2PSAction(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fFileName = "";
  fFile = 0;
  fViewport[0] = 0;
  fViewport[1] = 0;
  fViewport[2] = 0;
  fViewport[3] = 0;
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGL2PSAction::setLineWidth(
 int width
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  gl2psLineWidth( width );
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGL2PSAction::setPointSize(
 int size
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  gl2psPointSize( size );
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGL2PSAction::setViewport(
int a
,int b
,int winSizeX
,int winSizeY
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fViewport[0] = a;
  fViewport[1] = b;
  fViewport[2] = winSizeX;
  fViewport[3] = winSizeY;
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGL2PSAction::setFileName(
 const char* aFileName
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fFileName = aFileName;
}
//////////////////////////////////////////////////////////////////////////////
bool G4OpenGL2PSAction::enableFileWriting(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fFile = ::fopen(fFileName,"w");
  if(!fFile) {
    return false;
  }
  G4gl2psBegin();
  return true;
}
//////////////////////////////////////////////////////////////////////////////
void G4OpenGL2PSAction::disableFileWriting(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  gl2psEndPage();        
  ::fclose(fFile);
  fFile = 0;
}
//////////////////////////////////////////////////////////////////////////////
bool G4OpenGL2PSAction::fileWritingEnabled(
) const
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return (fFile?true:false);
}
//////////////////////////////////////////////////////////////////////////////
void G4OpenGL2PSAction::G4gl2psBegin(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fFile) return;
  int options = GL2PS_OCCLUSION_CULL | 
    GL2PS_BEST_ROOT | GL2PS_SILENT | GL2PS_DRAW_BACKGROUND | GL2PS_USE_CURRENT_VIEWPORT;
//   int options = GL2PS_OCCLUSION_CULL | 
//      GL2PS_BEST_ROOT | GL2PS_SILENT | GL2PS_DRAW_BACKGROUND;
  int sort = GL2PS_BSP_SORT;
  //int sort = GL2PS_SIMPLE_SORT;
  GLint buffsize = 0;
  buffsize += 1024*1024;
  
  glGetIntegerv(GL_VIEWPORT,fViewport);

  gl2psBeginPage("title","HEPVis::G4OpenGL2PSAction", 
                 fViewport,
                 GL2PS_EPS, 
                 sort, 
                 options, 
                 GL_RGBA,0, NULL,0,0,0,
                 buffsize, 
                 fFile,fFileName);    
}




#endif
