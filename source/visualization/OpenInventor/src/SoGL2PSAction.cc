// this :
#include "HEPVis/actions/SoGL2PSAction.h"

// Inventor :
#include <Inventor/elements/SoViewportRegionElement.h>
#include <Inventor/errors/SoDebugError.h>

#include "Geant4_gl2ps.h"

#include <stdio.h>

SO_ACTION_SOURCE(SoGL2PSAction)
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::initClass(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_ACTION_INIT_CLASS(SoGL2PSAction,SoGLRenderAction);
}
//////////////////////////////////////////////////////////////////////////////
SoGL2PSAction::SoGL2PSAction(
 const SbViewportRegion& aViewPortRegion
)
:SoGLRenderAction(aViewPortRegion)
,fFileName("out.ps")
,fEnable(FALSE)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_ACTION_CONSTRUCTOR(SoGL2PSAction);
}
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::setFileName(
 const char* aFileName
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!aFileName) return;
  fFileName = aFileName;
}
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::enableFileWriting(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fEnable = TRUE;
}
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::disableFileWriting(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fEnable = FALSE;
}
//////////////////////////////////////////////////////////////////////////////
SbBool SoGL2PSAction::fileWritingEnabled(
) const
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return fEnable;
}
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::beginTraversal(
 SoNode* aNode
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fEnable==FALSE) {

    SoGLRenderAction::beginTraversal(aNode);

  } else {
  
    const SbViewportRegion& vpr = getViewportRegion();

    SoViewportRegionElement::set(getState(),vpr);

    FILE* file = ::fopen(fFileName.getString(),"w");
    if(!file) {
      SoDebugError::post("SoGL2PSAction::beginTraversal",
                         "Cannot open file %s",fFileName.getString());
      return;
    }

    int options = GL2PS_OCCLUSION_CULL 
      | GL2PS_BEST_ROOT 
      | GL2PS_SILENT
      | GL2PS_DRAW_BACKGROUND;
    int sort = GL2PS_BSP_SORT;
    //int sort = GL2PS_SIMPLE_SORT;
    
    const SbVec2s& win = vpr.getWindowSize();
    GLint vp[4];
    vp[0] = 0;
    vp[1] = 0;
    vp[2] = win[0];
    vp[3] = win[1];

    int bufsize = 0;
    gl2psBeginPage("title","HEPVis::SoGL2PSAction", 
                   vp,
                   GL2PS_EPS, 
                   sort, 
                   options, 
                   GL_RGBA,0, NULL,0,0,0,
                   bufsize, 
                   file,fFileName.getString());
    
    traverse(aNode);
    
    /*int state = */
    gl2psEndPage();
        
    ::fclose(file);
  }
}
//////////////////////////////////////////////////////////////////////////////
SbBool SoGL2PSAction::addBitmap(
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
  if(fEnable==FALSE) return FALSE;
  GLboolean valid;
  glGetBooleanv(GL_CURRENT_RASTER_POSITION_VALID,&valid);
  if(valid==GL_FALSE) return FALSE;
  float pos[4];
  glGetFloatv(GL_CURRENT_RASTER_POSITION,pos);
  int xoff = -(int)(aXmove + aXorig);
  int yoff = -(int)(aYmove + aYorig);
  int x = (int)(pos[0] + xoff);
  int y = (int)(pos[1] + yoff);
  // Should clip against viewport area :
  GLint vp[4];
  glGetIntegerv(GL_VIEWPORT,vp);
  GLsizei w = aWidth;
  GLsizei h = aHeight;
  if(x+w>(vp[0]+vp[2])) w = vp[0]+vp[2]-x;
  if(y+h>(vp[1]+vp[3])) h = vp[1]+vp[3]-y;
  int s = 3 * w * h;
  if(s<=0) return FALSE;
  float* image = (float*)::malloc(s * sizeof(float));
  if(!image) return FALSE;
  glReadPixels(x,y,w,h,GL_RGB,GL_FLOAT,image);
  GLint status = gl2psDrawPixels(w,h,xoff,yoff,GL_RGB,GL_FLOAT,image);
  ::free(image);
  return (status!=GL2PS_SUCCESS ? FALSE : TRUE);
}
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::beginViewport(
)
/////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fEnable==FALSE) return;
  GLint vp[4];
  glGetIntegerv(GL_VIEWPORT,vp);
  gl2psBeginViewport(vp);
}
//////////////////////////////////////////////////////////////////////////////
void SoGL2PSAction::endViewport(
)
/////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fEnable==FALSE) return;
  gl2psEndViewport();
}
//////////////////////////////////////////////////////////////////////////////
SbBool SoGL2PSAction::addPixmap(
 float
,float
,int
,int
,float*
,SbBool
)
//////////////////////////////////////////////////////////////////////////////
// Dead born. Deprecated. Use addBitmap.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return FALSE;
}
