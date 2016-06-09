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
#ifdef G4VIS_BUILD_OI_DRIVER

/*----------------------------HEPVis----------------------------------------*/
/*                                                                          */
/* Node:             SoMarkerSet                                            */
/* Author:           Guy Barrand                                            */
/*                                                                          */
/*--------------------------------------------------------------------------*/

// this :
#include <HEPVis/nodes/SoMarkerSet.h>

#include <Inventor/errors/SoDebugError.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/nodes/SoPointSet.h>

#include <Inventor/elements/SoCoordinateElement.h>
#include <Inventor/elements/SoCacheElement.h>
#include <Inventor/elements/SoLazyElement.h>

#include <HEPVis/SbGL.h>
#include <HEPVis/actions/SoGL2PSAction.h>

static void drawMarker(SoAction*,int);
static GLubyte* getBitmap(int,int,char []); 

/* 
  "  x  "
  "  x  "
  "xxxxx"
  "  x  "
  "  x  "

  Should produce bitmap :
  0x20,0x20,0xf8,0x20,0x20

  The rows will be rendered down to top ; first row at bottom, last at top.
  In the below, '-' means that glBitmap will move the pointer to next byte.

  32103210 32103210 32103210 32103210 32103210
  ..1..--- ..1..--- 11111--- ..1..--- ..1..---

  0x20     0x20     0xf8     0x20     0x20
*/

///////////////////////////////////////////////////////////////
/// 5 5 ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
static char plus_5_5[] = {
  "  x  "
  "  x  "
  "xxxxx"
  "  x  "
  "  x  "
};
static char asterisk_5_5[] = {
  "x x x"
  " xxx "
  "  x  "
  " xxx "
  "x x x"
};
static char cross_5_5[] = {
  "x   x"
  " x x "
  "  x  "
  " x x "
  "x   x"
};
static char star_5_5[] = {
  "x x x"
  " xxx "
  "xxxxx"
  " xxx "
  "x x x"
};
static char circle_line_5_5[] = {
  " xxx "
  "x   x"
  "x   x"
  "x   x"
  " xxx "
};       
static char circle_filled_5_5[] = {
  " xxx "
  "xxxxx"
  "xxxxx"
  "xxxxx"
  " xxx "
};       
static char triangle_up_line_5_5[] = { //OpenGL will draw with y reversed.
  "xxxxx"
  " x x "
  " x x "
  "  x  "
  "  x  "
};
static char triangle_up_filled_5_5[] = {
  "xxxxx"
  " xxx "
  " xxx "
  "  x  "
  "  x  "
};
static char triangle_down_line_5_5[] = {
  "  x  "
  "  x  "
  " x x "
  " x x "
  "xxxxx"
};
static char triangle_down_filled_5_5[] = {
  "  x  "
  "  x  "
  " xxx "
  " xxx "
  "xxxxx"
};
static char david_star_line_5_5[] = {
  "  x  "
  "xxxxx"
  " x x "
  "xxxxx"
  "  x  "
};       
static char david_star_filled_5_5[] = {
  "  x  "
  "xxxxx"
  " xxx "
  "xxxxx"
  "  x  "
};       
static char swiss_cross_line_5_5[] = {
  " xxx "
  "xx xx"
  "x   x"
  "xx xx"
  " xxx "
};       
static char swiss_cross_filled_5_5[] = {
  " xxx "
  "xxxxx"
  "xxxxx"
  "xxxxx"
  " xxx "
};       
static char diamond_line_5_5[] = {
  "  x  "
  " x x "
  "x   x"
  " x x "
  "  x  "
};
static char diamond_filled_5_5[] = {
  "  x  "
  " xxx "
  "xxxxx"
  " xxx "
  "  x  "
};
static char square_line_5_5[] = {
  "xxxxx"
  "x   x"
  "x   x"
  "x   x"
  "xxxxx"
};
static char square_filled_5_5[] = {
  "xxxxx"
  "xxxxx"
  "xxxxx"
  "xxxxx"
  "xxxxx"
};
///////////////////////////////////////////////////////////////
/// 7 7 ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
static char plus_7_7[] = {
  "   x   "
  "   x   "
  "   x   "
  "xxxxxxx"
  "   x   "
  "   x   "
  "   x   "
};
static char asterisk_7_7[] = {
  "x  x  x"
  " x x x "
  "  xxx  "
  "   x   "
  "  xxx  "
  " x x x "
  "x  x  x"
};
static char cross_7_7[] = {
  "x     x"
  " x   x "
  "  xxx  "
  "   x   "
  "  xxx  "
  " x   x "
  "x     x"
};
static char star_7_7[] = {
  "x  x  x"
  " x x x "
  "  xxx  "
  "xxxxxxx"
  "  xxx  "
  " x x x "
  "x  x  x"
};
static char circle_line_7_7[] = {
  " xxxxx "
  "x     x"
  "x     x"
  "x     x"
  "x     x"
  "x     x"
  " xxxxx "
};       
static char circle_filled_7_7[] = {
  " xxxxx "
  "xxxxxxx"
  "xxxxxxx"
  "xxxxxxx"
  "xxxxxxx"
  "xxxxxxx"
  " xxxxx "
};       
static char triangle_up_line_7_7[] = { //OpenGL will draw with y reversed.
  "xxxxxxx"
  " x   x "
  " x   x "
  "  x x  "
  "  x x  "
  "   x   "
  "   x   "
};
static char triangle_up_filled_7_7[] = {
  "xxxxxxx"
  " xxxxx "
  " xxxxx "
  "  xxx  "
  "  xxx  "
  "   x   "
  "   x   "
};
static char triangle_down_line_7_7[] = {
  "   x   "
  "   x   "
  "  x x  "
  "  x x  "
  " x   x "
  " x   x "
  "xxxxxxx"
};
static char triangle_down_filled_7_7[] = {
  "   x   "
  "   x   "
  "  xxx  "
  "  xxx  "
  " xxxxx "
  " xxxxx "
  "xxxxxxx"
};
static char david_star_line_7_7[] = {
  "   x   "
  "xxxxxxx"
  " x   x "
  "  x x  "
  " x   x "
  "xxxxxxx"
  "   x   "
};       
static char david_star_filled_7_7[] = {
  "   x   "
  "xxxxxxx"
  " xxxxx "
  "  xxx  "
  " xxxxx "
  "xxxxxxx"
  "   x   "
};       
static char swiss_cross_line_7_7[] = {
  "  xxx  "
  "  x x  "
  "xxx xxx"
  "x     x"
  "xxx xxx"
  "  x x  "
  "  xxx  "
};       
static char swiss_cross_filled_7_7[] = {
  "  xxx  "
  "  xxx  "
  "xxxxxxx"
  "xxxxxxx"
  "xxxxxxx"
  "  xxx  "
  "  xxx  "
};       
static char diamond_line_7_7[] = {
  "   x   "
  "  x x  "
  " x   x "
  "x     x"
  " x   x "
  "  x x  "
  "   x   "
};
static char diamond_filled_7_7[] = {
  "   x   "
  "  xxx  "
  " xxxxx "
  "xxxxxxx"
  " xxxxx "
  "  xxx  "
  "   x   "
};
static char square_line_7_7[] = {
  "xxxxxxx"
  "x     x"
  "x     x"
  "x     x"
  "x     x"
  "x     x"
  "xxxxxxx"
};
static char square_filled_7_7[] = {
  "xxxxxxx"
  "xxxxxxx"
  "xxxxxxx"
  "xxxxxxx"
  "xxxxxxx"
  "xxxxxxx"
  "xxxxxxx"
};

///////////////////////////////////////////////////////////////
/// 9 9 ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
static char plus_9_9[] = {
  "    x    "
  "    x    "
  "    x    "
  "    x    "
  "xxxxxxxxx"
  "    x    "
  "    x    "
  "    x    "
  "    x    "
};
static char asterisk_9_9[] = {
  "x   x   x"
  " x  x  x "
  "  x x x  "
  "   xxx   "
  "    x    "
  "   xxx   "
  "  x x x  "
  " x  x  x "
  "x   x   x"
};
static char cross_9_9[] = {
  "x       x"
  " x     x "
  "  x   x  "
  "   x x   "
  "    x    "
  "   x x   "
  "  x   x  "
  " x     x "
  "x       x"
};
static char star_9_9[] = {
  "x   x   x"
  " x  x  x "
  "  x x x  "
  "   xxx   "
  "xxxxxxxxx"
  "   xxx   "
  "  x x x  "
  " x  x  x "
  "x   x   x"
};
static char circle_line_9_9[] = {
  "   xxx   "
  " xx   xx "
  " x     x "
  "x       x"
  "x       x"
  "x       x"
  " x     x "
  " xx   xx "
  "   xxx   "
};       
static char circle_filled_9_9[] = {
  "   xxx   "
  " xxxxxxx "
  " xxxxxxx "
  "xxxxxxxxx"
  "xxxxxxxxx"
  "xxxxxxxxx"
  " xxxxxxx "
  " xxxxxxx "
  "   xxx   "
};       
static char triangle_up_line_9_9[] = { //OpenGL will draw with y reversed.
  "xxxxxxxxx"
  " x     x "
  " x     x "
  "  x   x  "
  "  x   x  "
  "   x x   "
  "   x x   "
  "    x    "
  "    x    "
};
static char triangle_up_filled_9_9[] = {
  "xxxxxxxxx"
  " xxxxxxx "
  " xxxxxxx "
  "  xxxxx  "
  "  xxxxx  "
  "   xxx   "
  "   xxx   "
  "    x    "
  "    x    "
};
static char triangle_down_line_9_9[] = {
  "    x    "
  "    x    "
  "   x x   "
  "   x x   "
  "  x   x  "
  "  x   x  "
  " x     x "
  " x     x "
  "xxxxxxxxx"
};
static char triangle_down_filled_9_9[] = {
  "    x    "
  "    x    "
  "   xxx   "
  "   xxx   "
  "  xxxxx  "
  "  xxxxx  "
  " xxxxxxx "
  " xxxxxxx "
  "xxxxxxxxx"
};
static char david_star_line_9_9[] = {
  "    x    "
  "   x x   "
  "xxxxxxxxx"
  " x     x "
  "  x   x  "
  " x     x "
  "xxxxxxxxx"
  "   x x   "
  "    x    "
};       
static char david_star_filled_9_9[] = {
  "    x    "
  "   xxx   "
  "xxxxxxxxx"
  " xxxxxxx "
  "  xxxxx  "
  " xxxxxxx "
  "xxxxxxxxx"
  "   xxx   "
  "    x    "
};       
static char swiss_cross_line_9_9[] = {
  "   xxx   "
  "   x x   "
  "   x x   "
  "xxxx xxxx"
  "x       x"
  "xxxx xxxx"
  "   x x   "
  "   x x   "
  "   xxx   "
};       
static char swiss_cross_filled_9_9[] = {
  "   xxx   "
  "   xxx   "
  "   xxx   "
  "xxxxxxxxx"
  "xxxxxxxxx"
  "xxxxxxxxx"
  "   xxx   "
  "   xxx   "
  "   xxx   "
};       
static char diamond_line_9_9[] = {
  "    x    "
  "   x x   "
  "  x   x  "
  " x     x "
  "x       x"
  " x     x "
  "  x   x  "
  "   x x   "
  "    x    "
};
static char diamond_filled_9_9[] = {
  "    x    "
  "   xxx   "
  "  xxxxx  "
  " xxxxxxx "
  "xxxxxxxxx"
  " xxxxxxx "
  "  xxxxx  "
  "   xxx   "
  "    x    "
};
static char square_line_9_9[] = {
  "xxxxxxxxx"
  "x       x"
  "x       x"
  "x       x"
  "x       x"
  "x       x"
  "x       x"
  "x       x"
  "xxxxxxxxx"
};
static char square_filled_9_9[] = {
  "xxxxxxxxx"
  "xxxxxxxxx"
  "xxxxxxxxx"
  "xxxxxxxxx"
  "xxxxxxxxx"
  "xxxxxxxxx"
  "xxxxxxxxx"
  "xxxxxxxxx"
  "xxxxxxxxx"
};

static char* sFigures[54] = {
 plus_5_5,   //0
 asterisk_5_5,
 cross_5_5,
 star_5_5,
 circle_line_5_5,
 circle_filled_5_5,
 triangle_up_line_5_5,
 triangle_up_filled_5_5,
 triangle_down_line_5_5,
 triangle_down_filled_5_5,
 david_star_line_5_5,
 david_star_filled_5_5,
 swiss_cross_line_5_5,
 swiss_cross_filled_5_5,
 diamond_line_5_5,
 diamond_filled_5_5,
 square_line_5_5,
 square_filled_5_5, //17
 plus_7_7,
 asterisk_7_7,
 cross_7_7,
 star_7_7,
 circle_line_7_7,
 circle_filled_7_7,
 triangle_up_line_7_7,
 triangle_up_filled_7_7,
 triangle_down_line_7_7,
 triangle_down_filled_7_7,
 david_star_line_7_7,
 david_star_filled_7_7,
 swiss_cross_line_7_7,
 swiss_cross_filled_7_7,
 diamond_line_7_7,
 diamond_filled_7_7,
 square_line_7_7,
 square_filled_7_7, //35
 plus_9_9,
 asterisk_9_9,
 cross_9_9,
 star_9_9,
 circle_line_9_9,
 circle_filled_9_9,
 triangle_up_line_9_9,
 triangle_up_filled_9_9,
 triangle_down_line_9_9,
 triangle_down_filled_9_9,
 david_star_line_9_9,
 david_star_filled_9_9,
 swiss_cross_line_9_9,
 swiss_cross_filled_9_9,
 diamond_line_9_9,
 diamond_filled_9_9,
 square_line_9_9,
 square_filled_9_9 //53
};

SO_NODE_SOURCE(HEPVis_SoMarkerSet)
//////////////////////////////////////////////////////////////////////////////
void HEPVis_SoMarkerSet::initClass (
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_NODE_INIT_CLASS(HEPVis_SoMarkerSet,SoPointSet,"PointSet");
}
//////////////////////////////////////////////////////////////////////////////
HEPVis_SoMarkerSet::HEPVis_SoMarkerSet (
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_NODE_CONSTRUCTOR(HEPVis_SoMarkerSet);
  
  SO_NODE_ADD_FIELD(markerIndex,(CROSS_5_5));
}
//////////////////////////////////////////////////////////////////////////////
HEPVis_SoMarkerSet::~HEPVis_SoMarkerSet (
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}
//////////////////////////////////////////////////////////////////////////////
void HEPVis_SoMarkerSet::GLRender (
 SoGLRenderAction* aAction
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SoState* state = aAction->getState();

  const SoCoordinateElement* coordinateElement = 
    SoCoordinateElement::getInstance(state);
  if(coordinateElement==NULL) return;

  if(aAction->isOfType(SoGL2PSAction::getClassTypeId())) {
    SoCacheElement::invalidate(state);
  }

  const SbColor& color = SoLazyElement::getDiffuse(aAction->getState(),0);
  float red,green,blue;
  color.getValue(red,green,blue);

  int mark = markerIndex[0];

  int starti = startIndex.getValue();
  int pointn = numPoints.getValue();
  int pointi;

  glPushAttrib( (GLbitfield)(GL_CURRENT_BIT | GL_ENABLE_BIT));
  glDisable(GL_LIGHTING);
  glColor3f(red,green,blue);

#ifdef WIN32  
  //WIN32 : depth test is out over bitmap !
  glDisable(GL_DEPTH_TEST);
#endif

  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  for(pointi=starti;pointi<pointn;pointi++){
    const SbVec3f& vec = coordinateElement->get3(pointi);
    glRasterPos3f(vec[0],vec[1],vec[2]);
    // Do a push, pop to correct a deffect of Mesa-3.1. 
    // If not, further line drawing will have bad colors.
    // The glPopAttrib will compell a reinitialisation of
    // some internal Mesa state.
    //glPushAttrib(GL_ALL_ATTRIB_BITS);
    //glPopAttrib();
    //
    drawMarker(aAction,mark);
  }

  glPopAttrib();
}
//////////////////////////////////////////////////////////////////////////////
void drawMarker(
 SoAction* aAction
,int aStyle
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  GLsizei w = 0,h = 0;
  GLfloat xorig = 0,yorig = 0;
  GLfloat xmove = 0,ymove = 0;

  if((aStyle>=0)&&(aStyle<18)) {
    w = h = 5;
    xorig = yorig = 2;
    GLubyte* bitmap = getBitmap(w,h,sFigures[aStyle]);
    glBitmap(w,h,xorig,yorig,0.,0.,bitmap);
    delete [] bitmap;
  } else if((aStyle>=18)&&(aStyle<36)) {
    w = h = 7;
    xorig = yorig = 3;
    GLubyte* bitmap = getBitmap(w,h,sFigures[aStyle]);
    glBitmap(w,h,xorig,yorig,0.,0.,bitmap);
    delete [] bitmap;
  } else if((aStyle>=36)&&(aStyle<54)) {
    w = h = 9;
    xorig = yorig = 4;
    GLubyte* bitmap = getBitmap(w,h,sFigures[aStyle]);
    glBitmap(w,h,xorig,yorig,0.,0.,bitmap);
    delete [] bitmap;
  } else {
    return;
  }

  if(aAction->isOfType(SoGL2PSAction::getClassTypeId())) {
    ((SoGL2PSAction*)aAction)->addBitmap(w,h,xorig,yorig,xmove,ymove);
  }

}
//////////////////////////////////////////////////////////////////////////////
GLubyte* getBitmap(
 int aW
,int aH
,char aFigure[]
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  int index = 0;
  GLubyte* bitmap = new GLubyte[aW * aH + 1];
  int ichar = 0;
  int ibit = 0;
  unsigned char byte = 0;
  for ( int row = 0; row < aH; row++ ){
    for ( int col = 0; col < aW; col++){ 
      unsigned char c = aFigure[ichar];
      ichar++;
      if(c==' ') {
        ibit++;
      } else { 
        byte += (1<<(7-ibit));
        ibit++;
      }
      if(ibit==8) {
        //unsigned char h = byte / 16;
        //unsigned char l = byte % 16;
        //printf("0x%x%x\n",h,l);
        bitmap[index] = byte;
        index++;
        ibit = 0;
        byte = 0;
      }

    }
    if(ibit!=8) { //Jump to next byte.
      //unsigned char h = byte / 16;
      //unsigned char l = byte % 16;
      //printf("0x%x%x\n",h,l);
      bitmap[index] = byte;
      index++;
      ibit = 0;
      byte = 0;
    }
  }
  return bitmap; 
}

#endif
