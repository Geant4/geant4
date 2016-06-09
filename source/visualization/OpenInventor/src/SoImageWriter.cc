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
/* Node:             SoImageWriter                                          */
/* Author:           Guy Barrand                                            */
/*                                                                          */
/*--------------------------------------------------------------------------*/

// this :
#include <HEPVis/nodes/SoImageWriter.h>

#include <Inventor/errors/SoDebugError.h>
#include <Inventor/elements/SoViewportRegionElement.h>
#include <Inventor/actions/SoGLRenderAction.h>

#include <HEPVis/SbGL.h>
#include <HEPVis/SbPainterPS.h>
//#include <HEPVis/SbGIF.h>

#include <stdlib.h>

typedef struct {
  unsigned char red;
  unsigned char green;
  unsigned char blue;
} Pixel;
typedef unsigned char Uchar;

//static void getImagePixels(int,int,float*,int&,
//                         Uchar*&,Uchar*&,Uchar*&,Uchar*&);

static int sWidth = 0;
static int sHeight = 0;
static float* sImage = 0;
static int getRGB(unsigned int,unsigned int,double&,double&,double&);

SO_NODE_SOURCE(SoImageWriter)
//////////////////////////////////////////////////////////////////////////////
void SoImageWriter::initClass (
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_NODE_INIT_CLASS(SoImageWriter,SoNode,"Node");
}
//////////////////////////////////////////////////////////////////////////////
SoImageWriter::SoImageWriter(
)
:fEnabled(FALSE)
,fStatus(FALSE)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_NODE_CONSTRUCTOR(SoImageWriter);
  //SO_NODE_ADD_FIELD(format,(POST_SCRIPT));
  SO_NODE_ADD_FIELD(fileName,("out.ps"));
  
  //SO_NODE_DEFINE_ENUM_VALUE(Format,POST_SCRIPT);
  //SO_NODE_DEFINE_ENUM_VALUE(Format,GIF);
  
  //SO_NODE_SET_SF_ENUM_TYPE(format,Format);
}
//////////////////////////////////////////////////////////////////////////////
SoImageWriter::~SoImageWriter (
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}
//////////////////////////////////////////////////////////////////////////////
void SoImageWriter::enable(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fEnabled = TRUE;
}
//////////////////////////////////////////////////////////////////////////////
void SoImageWriter::disable(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fEnabled = FALSE;
}
//////////////////////////////////////////////////////////////////////////////
SbBool SoImageWriter::getStatus(
) const
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return fStatus;
}
//////////////////////////////////////////////////////////////////////////////
void SoImageWriter::GLRender(
 SoGLRenderAction* aAction
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fStatus = FALSE;
  //printf("debug : SoImageWriter::GLRender : enabled : %d\n",fEnabled);
  if(!fEnabled) return;
  SbViewportRegion vpr = SoViewportRegionElement::get(aAction->getState());
  const SbVec2s& win = vpr.getWindowSize();
  int w = win[0];
  int h = win[1];
  if((w*h)<=0) {
    SoDebugError::postInfo("SoImageWriter::GLRender","null area window !");
    return;
  }

  int x = 0;
  int y = 0;
  int s = 3 * w * h;
  float* image = new float[s];
  if(!image) return;

  //printf("debug : SoImageWriter::GLRender : %d %d %d %d\n",x,y,w,h);

  //glReadPixels(x,y,w,h,GL_RGB,GL_UNSIGNED_BYTE,image); Don't work !
  glReadPixels(x,y,w,h,GL_RGB,GL_FLOAT,image);

  //Format fm = (Format)format.getValue();
  //if(fm==GIF) {  
/*
    FILE* file = fopen(fileName.getValue().getString(),"wb");
    if(!file) {
      SoDebugError::postInfo("SoImageWriter::GLRender",
        "can't open file \"%s\".",fileName.getValue().getString());
    } else {
      int coln;
      Uchar* rs;
      Uchar* gs;
      Uchar* bs;
      Uchar* data;
      getImagePixels(w,h,image,coln,rs,gs,bs,data);
      
      SbGIF::putBytesInStream(file,data,w,h,coln,rs,gs,bs);
      
      delete [] data;

      if(rs) free(rs);
      if(gs) free(gs);
      if(bs) free(bs);

      fclose(file);

      fStatus = TRUE;
    }
  } else {
*/

    SbPainterPS painterPS;
    painterPS.openFileForWriting(fileName.getValue().getString());
    if(!painterPS.getStream()) {
      SoDebugError::postInfo("SoImageWriter::GLRender",
        "can't open file \"%s\".",fileName.getValue().getString());
    } else {
      painterPS.setWindowSize(w,h);
      //painterPS.setBitsPerPixel(8);
      painterPS.setBitsPerPixel(4);
      painterPS.beginTraversal();
      painterPS.clearColorBuffer(1.,1.,1.);

      sWidth = w;
      sHeight = h;
      sImage = image;
      painterPS.putImageInStream((unsigned int)w,(unsigned int)h,getRGB);
      
      painterPS.endTraversal();
      
      painterPS.closeStream();
      
      fStatus = TRUE;
    }
  //}
  delete [] image;

}
/*
//////////////////////////////////////////////////////////////////////////////
void getImagePixels (
 int aWidth
,int aHeight
,float* aImage
,int& aColorn
,Uchar*& aReds
,Uchar*& aGreens
,Uchar*& aBlues
,Uchar*& aData
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  aColorn = 0;
  aReds = 0;
  aGreens = 0;
  aBlues = 0;
  aData = 0;
  if( (aWidth * aHeight) <=0) return;
  int size = 256;
  Uchar* rs = (Uchar*)malloc(size * sizeof(Uchar));
  Uchar* gs = (Uchar*)malloc(size * sizeof(Uchar));
  Uchar* bs = (Uchar*)malloc(size * sizeof(Uchar));
  Uchar* data = new Uchar[aWidth * aHeight];
  if( !rs || !gs || !bs || !data ) {
    if(rs) free(rs);
    if(gs) free(gs);
    if(bs) free(bs);
    delete [] data;
    return;
  }
  int pixeln = 0;
  int row,col;
  Uchar red,green,blue;
  Uchar ored = 0,ogreen = 0,oblue = 0;
  float* pimag = aImage;
  Uchar* pdata = 0;
  Uchar index = 0;
  int status = 0;
  for(row=0;row<aHeight;row++) {
    pdata = data + (aHeight - 1 - row) * aWidth;
    for(col=0;col<aWidth;col++){ 
      red   = (Uchar)(255 * (*pimag));pimag++;
      green = (Uchar)(255 * (*pimag));pimag++;
      blue  = (Uchar)(255 * (*pimag));pimag++;
      //printf("debug : %d %d : %d %d %d\n",row,col,red,green,blue);
      if( (pixeln==0) || (red!=ored) || (green!=ogreen) || (blue!=oblue) ){
        // Search exact color :
        int found = 0;
        for(int count=0;count<pixeln;count++){
          if( (red==rs[count]) && (green==gs[count]) && (blue==bs[count]) ){
            found = 1;
            index = count;
            break;
          }
        }
        if(found==0){
          if(pixeln>=256) {
            // We can't store more than 256 on an Uchar.
            // Search closest color :
            int dr,dg,db;
            int PRECISION = 20;
            int closest = 0;
            for(int count=0;count<pixeln;count++){
              dr = red   - rs[count];dr = dr<0 ? -dr : dr;
              dg = green - gs[count];dg = dg<0 ? -dg : dg;
              db = blue  - bs[count];db = db<0 ? -db : db;
              if( (dr<=PRECISION) && (dg<=PRECISION) && (db<=PRECISION) ){
                closest = 1;
                index = count;
                break;
              }
            }
            if(closest==0) {
              index = 0;
              status  = 1;
            }
          } else {
            if(pixeln>=size){
              size += 256;
              rs = (Uchar*)realloc(rs,size * sizeof(Uchar));
              gs = (Uchar*)realloc(gs,size * sizeof(Uchar));
              bs = (Uchar*)realloc(bs,size * sizeof(Uchar));
              if( !rs || !gs || !bs ) {
                if(rs) free(rs);
                if(gs) free(gs);
                if(bs) free(bs);
                delete [] data;
                return;
              } 
            } 
            //printf("debug : SoImageWriter pixeln %d : %d %d %d\n",
            //   pixeln,red,green,blue);
            rs[pixeln] = red;
            gs[pixeln] = green;
            bs[pixeln] = blue;
            index = pixeln;
            pixeln++;
          }
        }
      }
      *pdata = index;
      pdata++;
      ored = red;
      ogreen = green;
      oblue = blue;
    }
  }
  if(status==1) 
    printf("SoImageWriter : more than 256 colors in picture ; some colors approximated.\n");
  aColorn = pixeln;
  aReds = rs;
  aGreens = gs;
  aBlues = bs;
  aData = data;
}
*/
//////////////////////////////////////////////////////////////////////////////
int getRGB(
 unsigned int aX
,unsigned int aY
,double& aRed
,double& aGreen
,double& aBlue
)
//////////////////////////////////////////////////////////////////////////////
// OpenGL image is from down to up.
// PS image is up to down.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  float* pimag = sImage + 3 * (sWidth * (sHeight - 1 - aY) + aX);
  aRed   = *pimag;pimag++;
  aGreen = *pimag;pimag++;
  aBlue  = *pimag;pimag++;
  return 1;
}

#endif
