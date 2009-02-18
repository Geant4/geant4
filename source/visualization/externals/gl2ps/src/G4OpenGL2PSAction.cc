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
// $Id: G4OpenGL2PSAction.cc,v 1.1 2009-02-18 09:54:12 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#define G4DEBUG_VIS_OGL 1

#ifdef G4VIS_BUILD_OPENGL_DRIVER
 #define G4VIS_BUILD_OPENGL_GL2PS 
#endif
#ifdef G4VIS_BUILD_OI_DRIVER
 #define G4VIS_BUILD_OPENGL_GL2PS 
#endif

#ifdef G4VIS_BUILD_OPENGL_GL2PS

#include "G4OpenGL2PSAction.hh"
#include "Geant4_gl2ps.h"


G4OpenGL2PSAction::G4OpenGL2PSAction(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGL2PSAction::G4OpenGL2PSAction\n");
#endif
  fFileName = "";
  fFile = 0;
  fViewport[0] = 0;
  fViewport[1] = 0;
  fViewport[2] = 0;
  fViewport[3] = 0;
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGL2PSAction::setFileName(
 const char* aFileName
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGL2PSAction::setFileName\n");
#endif
  fFileName = aFileName;
}
//////////////////////////////////////////////////////////////////////////////
bool G4OpenGL2PSAction::enableFileWriting(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGL2PSAction::enabledFileWriting\n");
#endif
  fFile = ::fopen(fFileName,"w");
  if(!fFile) {
    return false;
  }
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGL2PSAction::enabledFileWriting-------------\n");
#endif
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
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGL2PSAction::fileWrintingEnabled\n");
#endif
  return (fFile?true:false);
}
//////////////////////////////////////////////////////////////////////////////
void G4OpenGL2PSAction::G4gl2psBegin(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGL2PSAction::G4gl2psBegin\n");
#endif
  if(!fFile) return;
  int options = GL2PS_OCCLUSION_CULL | 
     GL2PS_BEST_ROOT | GL2PS_SILENT | GL2PS_DRAW_BACKGROUND;
  int sort = GL2PS_BSP_SORT;
  //int sort = GL2PS_SIMPLE_SORT;
    
  glGetIntegerv(GL_VIEWPORT,fViewport);

  int bufsize = 0;
  gl2psBeginPage("title","HEPVis::G4OpenGL2PSAction", 
                 fViewport,
                 GL2PS_EPS, 
                 sort, 
                 options, 
                 GL_RGBA,0, NULL,0,0,0,
                 bufsize, 
                 fFile,fFileName);    
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGL2PSAction::G4gl2psBegin END\n");
#endif
}

#endif
