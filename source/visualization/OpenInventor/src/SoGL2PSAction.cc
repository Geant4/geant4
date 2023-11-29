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

/*----------------------------HEPVis----------------------------------------*/
/*                                                                          */
/* Node:             SoGL2PSAction                                          */
/* Author:           Guy Barrand                                            */
/*                                                                          */
/*--------------------------------------------------------------------------*/

// this :
#include <HEPVis/actions/SoGL2PSAction.h>

// Inventor :
#include <Inventor/elements/SoViewportRegionElement.h>
#include <Inventor/errors/SoDebugError.h>

#include <Inventor/system/gl.h>

#include <stdio.h>

SO_ACTION_SOURCE(SoGL2PSAction)
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::initClass(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  static bool first = true;
  if (first) {
    first = false;
    SO_ACTION_INIT_CLASS(SoGL2PSAction,SoGLRenderAction);
  }
}
//////////////////////////////////////////////////////////////////////////////
SoGL2PSAction::SoGL2PSAction(
 const SbViewportRegion& aViewPortRegion
)
:SoGLRenderAction(aViewPortRegion)
,fContext(0)
,fFile(0)
,fFileName("out.pdf")
,fTitle("title")
,fProducer("HEPVis::SoGL2PSAction")
,fFormat(TOOLS_GL2PS_PDF)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_ACTION_CONSTRUCTOR(SoGL2PSAction);
}

//////////////////////////////////////////////////////////////////////////////
SoGL2PSAction::~SoGL2PSAction()
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  closeFile();
}

//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::setFileName(const std::string& aFileName)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fFileName = aFileName;
}

//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::setTitleAndProducer(const std::string& aTitle,const std::string& aProducer)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fTitle = aTitle;
  fProducer = aProducer;
}

void SoGL2PSAction::setExportImageFormat_PS()  {fFormat = TOOLS_GL2PS_PS;}
void SoGL2PSAction::setExportImageFormat_EPS() {fFormat = TOOLS_GL2PS_EPS;}
void SoGL2PSAction::setExportImageFormat_TEX() {fFormat = TOOLS_GL2PS_TEX;}
void SoGL2PSAction::setExportImageFormat_PDF() {fFormat = TOOLS_GL2PS_PDF;}
void SoGL2PSAction::setExportImageFormat_SVG() {fFormat = TOOLS_GL2PS_SVG;}
void SoGL2PSAction::setExportImageFormat_PGF() {fFormat = TOOLS_GL2PS_PGF;}

//////////////////////////////////////////////////////////////////////////////
bool SoGL2PSAction::enableFileWriting()
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!openFile()) {
    SoDebugError::post("SoGL2PSAction::enableFileWriting",
                       "openFile() failed for fil %s",
                       fFileName.c_str());
    return false;
  }
#ifdef __COIN__
#else //SGI
  const SbViewportRegion& vpr = getViewportRegion();
  SoViewportRegionElement::set(getState(),vpr);
  SbVec2s origin = vpr.getViewportOriginPixels();
  SbVec2s size = vpr.getViewportSizePixels();
  if(!beginPage(origin[0],origin[1],size[0],size[1])) {
    SoDebugError::post("SoGL2PSAction::enableFileWriting","beginPage() failed");
    return false;
  }
#endif
  return true;
}
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::disableFileWriting(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef __COIN__
#else //SGI
  endPage();
#endif
  closeFile();
}

//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::beginTraversal(
 SoNode* aNode
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fContext && fFile) {
#ifdef __COIN__
    const SbViewportRegion& vpr = getViewportRegion();
    SoViewportRegionElement::set(getState(),vpr);
    SbVec2s origin = vpr.getViewportOriginPixels();
    SbVec2s size = vpr.getViewportSizePixels();
    if(!beginPage(origin[0],origin[1],size[0],size[1])) {
      SoDebugError::post("SoGL2PSAction::beginTraversal","beginPage() failed");
      return;
    }
    traverse(aNode);
    if(!endPage()) {
      SoDebugError::post("SoGL2PSAction::beginTraversal","endPage() failed");
      return;
    }
#else //SGI
    SoGLRenderAction::beginTraversal(aNode);
#endif
  } else {
    SoGLRenderAction::beginTraversal(aNode);
  }
}

#include <tools/gl2ps>

//////////////////////////////////////////////////////////////////////////////
bool SoGL2PSAction::openFile()
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile) {
    ::fclose(fFile);
    fFile = 0;
  }
  if(fContext) {
    ::tools_gl2psDeleteContext(fContext);
    fContext = 0;
  }

  fContext = ::tools_gl2psCreateContext();
  if(!fContext) return false;

  tools_gl2ps_gl_funcs_t _funcs = {
    (tools_glIsEnabled_func)glIsEnabled,
    (tools_glBegin_func)glBegin,
    (tools_glEnd_func)glEnd,
    (tools_glGetFloatv_func)glGetFloatv,
    (tools_glVertex3f_func)glVertex3f,
    (tools_glGetBooleanv_func)glGetBooleanv,
    (tools_glGetIntegerv_func)glGetIntegerv,
    (tools_glRenderMode_func)glRenderMode,
    (tools_glFeedbackBuffer_func)glFeedbackBuffer,
    (tools_glPassThrough_func)glPassThrough
  };
  ::tools_gl2ps_set_gl_funcs(fContext,&_funcs);
  
  fFile = ::fopen(fFileName.c_str(),"wb");
  if(!fFile) {
    ::tools_gl2psDeleteContext(fContext);
    fContext = 0;
    return false;
  }

  return true;
}

//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::closeFile()
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile) {
    ::fclose(fFile);
    fFile = 0;
  }
  if(fContext) {
    ::tools_gl2psDeleteContext(fContext);
    fContext = 0;
  }    
}

//////////////////////////////////////////////////////////////////////////////
bool SoGL2PSAction::beginPage(int a_x,int a_y,int a_w,int a_h)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fContext) return false;
  if(!fFile) return false;

  if( (a_w<=0) || (a_h<=0) ) return false;

  int options = 
    TOOLS_GL2PS_BEST_ROOT |
    TOOLS_GL2PS_DRAW_BACKGROUND |
    TOOLS_GL2PS_USE_CURRENT_VIEWPORT;
  int sort = TOOLS_GL2PS_BSP_SORT;

  int vp[4];
  vp[0] = a_x;
  vp[1] = a_y;
  vp[2] = a_w;
  vp[3] = a_h;

  int bufferSize = 0;
  
  tools_GLint res = ::tools_gl2psBeginPage
    (fContext,fTitle.c_str(),fProducer.c_str(),
     vp,fFormat,sort,options,TOOLS_GL_RGBA,0, NULL,0,0,0,
     bufferSize,fFile,fFileName.c_str());
  if (res == TOOLS_GL2PS_ERROR) return false;

  // enable blending for all
  ::tools_gl2psEnable(fContext,TOOLS_GL2PS_BLEND);
  
  return true;
}

//////////////////////////////////////////////////////////////////////////////
bool SoGL2PSAction::endPage()
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  int _status = 0;
  if(fContext) {
    _status = ::tools_gl2psEndPage(fContext);
  }
  if (_status == TOOLS_GL2PS_OVERFLOW) return false;
  return true;
}

//////////////////////////////////////////////////////////////////////////////
bool SoGL2PSAction::addBitmap(
 int aWidth
,int aHeight
,float aXorig
,float aYorig
,float aXmove
,float aYmove
)
/////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fContext) return false;
  GLboolean valid;
  ::glGetBooleanv(GL_CURRENT_RASTER_POSITION_VALID,&valid);
  if(!valid) return false;
  float pos[4];
  ::glGetFloatv(GL_CURRENT_RASTER_POSITION,pos);
  int xoff = -(int)(aXmove + aXorig);
  int yoff = -(int)(aYmove + aYorig);
  int x = (int)(pos[0] + xoff);
  int y = (int)(pos[1] + yoff);
  // Should clip against viewport area :
  GLint vp[4];
  ::glGetIntegerv(GL_VIEWPORT,vp);
  GLsizei w = aWidth;
  GLsizei h = aHeight;
  if(x+w>(vp[0]+vp[2])) w = vp[0]+vp[2]-x;
  if(y+h>(vp[1]+vp[3])) h = vp[1]+vp[3]-y;
  int s = 3 * w * h;
  if(s<=0) return false;
  float* image = (float*)::malloc(s * sizeof(float));
  if(!image) return false;
  ::glReadPixels(x,y,w,h,GL_RGB,GL_FLOAT,image);
  GLint status = ::tools_gl2psDrawPixels(fContext,w,h,xoff,yoff,GL_RGB,GL_FLOAT,image);
  ::free(image);
  return (status!=TOOLS_GL2PS_SUCCESS ? false : true);
}

