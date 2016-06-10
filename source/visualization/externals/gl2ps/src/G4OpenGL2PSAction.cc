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
// $Id: G4OpenGL2PSAction.cc 88206 2015-02-03 13:01:27Z gcosmo $
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

#include <limits>
#include <cstdlib>
#include <cstring>

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
  fBufferSize = 2048;
  fBufferSizeLimit = (std::numeric_limits<GLint>::max)();
  fExportImageFormat = GL2PS_PDF;
  resetBufferSizeParameters();
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGL2PSAction::resetBufferSizeParameters(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fBufferSize = 2048;
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
    fFileName = (strncpy((char*)malloc((unsigned)strlen(aFileName) + 1), aFileName, (unsigned)strlen(aFileName) + 1));
}
//////////////////////////////////////////////////////////////////////////////
bool G4OpenGL2PSAction::enableFileWriting(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fFile = ::fopen(fFileName,"wb");
  if(!fFile) {
    return false;
  }

  // No buffering for output file
  setvbuf ( fFile , NULL , _IONBF , 2048 );
  return G4gl2psBegin();
}
//////////////////////////////////////////////////////////////////////////////
bool G4OpenGL2PSAction::disableFileWriting(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  int state = gl2psEndPage();
  ::fclose(fFile);
  if (state == GL2PS_OVERFLOW) {
    return false;
  }
  fFile = 0;
  return true;
}
//////////////////////////////////////////////////////////////////////////////
bool G4OpenGL2PSAction::extendBufferSize(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  // extend buffer size *2
  if (fBufferSize < (fBufferSizeLimit/2)) {
    fBufferSize = fBufferSize*2;
    return true;
  }
  return false;
}


// FWJ
void G4OpenGL2PSAction::setBufferSize(int newSize)
{
  fBufferSize = (newSize < int(fBufferSizeLimit))
              ? GLint(newSize) : fBufferSizeLimit;
}


bool G4OpenGL2PSAction::fileWritingEnabled(
) const
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return (fFile?true:false);
}
//////////////////////////////////////////////////////////////////////////////
bool G4OpenGL2PSAction::G4gl2psBegin(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fFile) return false;
  int options = 
    GL2PS_BEST_ROOT | GL2PS_DRAW_BACKGROUND |GL2PS_USE_CURRENT_VIEWPORT;
  int sort = GL2PS_BSP_SORT;

  glGetIntegerv(GL_VIEWPORT,fViewport);
  GLint res = gl2psBeginPage("Geant4 output","Geant4",
                 fViewport,
                 fExportImageFormat,
                 sort, 
                 options, 
                 GL_RGBA,0, NULL,0,0,0,
                 fBufferSize,
                 fFile,fFileName);    
  if (res == GL2PS_ERROR) {
    return false;
  }
  // enable blending for all
  gl2psEnable(GL2PS_BLEND);
  
  return true;
}

void G4OpenGL2PSAction::setExportImageFormat(unsigned int type){
  if(!fFile) {
    fExportImageFormat = type;
  } else {
    //  Could not change the file type at this step. Please change it before enableFileWriting()
  }
}




#endif
