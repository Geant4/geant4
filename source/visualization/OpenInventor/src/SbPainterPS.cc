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
/* Node:             SbPainterPS                                            */
/* Author:           Guy Barrand                                            */
/*                                                                          */
/*--------------------------------------------------------------------------*/
// this :
#include <HEPVis/SbPainterPS.h>

//#include <HEPVis/SbString.h>
#define STRDUP(str)  ((str) != NULL ? (::strcpy((char*)::malloc((unsigned)::strlen(str) + 1), str)) : (char*)NULL)
#define STRDEL(str) {if((str)!=NULL) {::free(str);str=NULL;}}

//#define DEBUG
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <locale.h>

#define METAFILE_DEFAULT "out.ps"
#define METAFILE_SCALE 1.

static char* GetDate();
static double ConvertRGB_ToGrey(double,double,double);
//////////////////////////////////////////////////////////////////////////////
SbPainterPS::SbPainterPS(
)
:fDeviceWidth((8.5-1.) * 72. * METAFILE_SCALE) /* 540. * METAFILE_SCALE */
,fDeviceHeight(11.     * 72. * METAFILE_SCALE) /* 792. * METAFILE_SCALE */
,fPageNumber(0)
,fMarkerSize(2.)
,fFile(NULL)
,fFileName(NULL)
,fGSave(0)
,fBufferCount(0)
,fBufferString(NULL)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fParams.shade = Color;
  fParams.portrait = 1;
  fParams.nbit = 2;
  fParams.doBack = 1;
  fParams.lineWidth = -1.;
  fBufferPointer[0] = '\0';
#ifdef WIN32
  ::setlocale(LC_NUMERIC,"USA");
#endif
}
//////////////////////////////////////////////////////////////////////////////
SbPainterPS::~SbPainterPS(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile!=NULL) closeStream ();
  if(fBufferString!=NULL) ::free(fBufferString);
  fBufferString = NULL;
  if(fGSave!=0) {
    ::printf("SbPainterPS : bad gsave/grestore balance : %d.\n",fGSave);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::beginTraversal (
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile==NULL) openFileForWriting(NULL);
  if(fFile==NULL) return;
  putBeginPageInStream();
  putPageScaleInStream((float)fWindowWidth,(float)fWindowHeight);
  putSaveStateInStream();
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::endTraversal(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile==NULL) return;
  putFrameInStream(0.0,0.0,0.0,(float)fWindowWidth,(float)fWindowHeight);
  putRestoreStateInStream();
  putEndPageInStream();
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::clearColorBuffer(
 float aRed
,float aGreen
,float aBlue
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile==NULL) return;
  putBackgroundInStream(aRed,aGreen,aBlue,
                           (float)fWindowWidth,(float)fWindowHeight);
}
/*
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::drawPrimitive (
 SbPrimitiveType aType
,int aPointn
,float* aXs
,float* aYs
,float* //aZs
,const SbPainterContext& aAtb
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile==NULL) return;
  switch(aType) {
  case SbPrimitivePoints:
    drawMarkers(aPointn,
                aXs,aYs,
                aAtb.fRed,aAtb.fGreen,aAtb.fBlue,
                aAtb.fMarkerStyle,aAtb.fMarkerSize);
    break;
  case SbPrimitiveLineStrip:
  case SbPrimitiveLineLoop:
    drawLines(aPointn,
              aXs,aYs,
              aAtb.fRed,aAtb.fGreen,aAtb.fBlue,
              aAtb.fLineStyle,aAtb.fLineWidth);
    break;
  case SbPrimitivePolygon:
    drawPolygon(aPointn,
                aXs,aYs,
                aAtb.fRed,aAtb.fGreen,aAtb.fBlue,
                aAtb.fAreaStyle);
    break;
  default:
    break;
  }
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::drawPolygon(
 int aPointn
,float* aXs
,float* aYs
,float aRed
,float aGreen
,float aBlue
,const SbAreaStyle& //aStyle
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile==NULL) return;
  if(aPointn<=0) return;
  putNewPathInStream(); 
  putMoveInStream(aXs[0],aYs[0]);
  for(int count=1;count<aPointn;count++) {
    putLineToInStream(aXs[count] - aXs[count-1],
                      aYs[count] - aYs[count-1]);
  }
  if ( (aXs[0]==aXs[aPointn-1]) &&
       (aYs[0]==aYs[aPointn-1]) ) 
    putClosePathInStream();
  putRGB_InStream(aRed,aGreen,aBlue);
  putFillInStream();
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::drawLines(
 int aPointn
,float* aXs
,float* aYs
,float aRed
,float aGreen
,float aBlue
,const SbLineStyle& aStyle
,int aWidth
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile==NULL) return;
  if(aPointn<=0) return;
  putMoveInStream(aXs[0],aYs[0]);
  for(int count=1;count<aPointn;count++) {
    putLineToInStream(aXs[count] - aXs[count-1],
                      aYs[count] - aYs[count-1]);
  }
  if ( (aXs[0]==aXs[aPointn-1]) &&
       (aYs[0]==aYs[aPointn-1]) ) 
    putClosePathInStream();
  putRGB_InStream(aRed,aGreen,aBlue);
  putLineWidthInStream(aWidth);
  putCapInStream(1);
  putLineStyleInStream(aStyle);
  putStrokeInStream();
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::drawMarkers (
 int aPointn
,float* aXs
,float* aYs
,float aRed
,float aGreen
,float aBlue
,const SbMarkerStyle& aStyle 
,int aSize
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile==NULL) return;
  float mark_size  = (float)(aSize <=0 ? 1. : aSize);   
  mark_size *= 0.6F;
  if(aStyle==SbMarkerCircleLine) {
    putNewPathInStream();
    int icount = 1;
    for(int count=0;count<aPointn;count++) {
      putCircleInStream(aXs[count],aYs[count],mark_size);
#define MAX_PATH_POINT 100
      if(icount==MAX_PATH_POINT) {
        putRGB_InStream(aRed,aGreen,aBlue);
        putLineWidthInStream(1);
        putCapInStream(1);
        putStrokeInStream();
        icount = 1;
        if(count!=aPointn-1) putNewPathInStream();
      } else {
        icount++;
      }
    }
    putRGB_InStream(aRed,aGreen,aBlue);
    putLineWidthInStream(1);
    putCapInStream(1);
    putStrokeInStream();
  } else {
    putNewPathInStream();
    int icount = 1;
    for(int count=0;count<aPointn;count++) {
      putMoveInStream(aXs[count],aYs[count]);
      putMarkerSizeInStream(mark_size);
      putMarkerStyleInStream(aStyle);
      if(icount==MAX_PATH_POINT) {
        putRGB_InStream(aRed,aGreen,aBlue);
        putLineWidthInStream(1);
        putCapInStream(1);
        putStrokeInStream();
        icount = 1;
        if(count!=aPointn-1) putNewPathInStream();
      } else {
        icount++;
      }
    }
    putRGB_InStream(aRed,aGreen,aBlue);
    putLineWidthInStream(1);
    putCapInStream(1);
    putStrokeInStream();
  }
}
*/
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::setColorScheme(
 int aShade
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fParams.shade = aShade;
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::setOrientation(
 int aPortrait
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fParams.portrait = aPortrait;
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::setBackgroundDrawn(
 int aDoback
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fParams.doBack = aDoback;
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::setBitsPerPixel(
 int aNbit
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if( (aNbit==2) || (aNbit==4) || (aNbit==8) )
    fParams.nbit = aNbit;
  else 
    fParams.nbit = 2;
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::setLineWidth(
 int aWidth
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fParams.lineWidth = (float)aWidth;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::setFileName(
 const char* aString
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  STRDEL(fFileName);
  fFileName = STRDUP(aString);
}
//////////////////////////////////////////////////////////////////////////////
const char* SbPainterPS::getFileName(
) const
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return fFileName;
}
//////////////////////////////////////////////////////////////////////////////
void* SbPainterPS::getStream(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return fFile;
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::openFileForWriting(
 const char* aString
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile!=NULL) closeStream ();
  if( (aString==NULL) || (*aString=='\0') ) {
    if( (fFileName==NULL) || (*fFileName=='\0') ) { // Take default name :
      fFile = ::fopen(METAFILE_DEFAULT,"wb");
      STRDEL(fFileName);
      fFileName = STRDUP(METAFILE_DEFAULT);
    } else {
      fFile = ::fopen(fFileName,"wb");
    }
  } else {
    fFile = ::fopen(aString,"wb");
    STRDEL(fFileName);
    fFileName = STRDUP(aString);
  }
  if(fFile==NULL) return;

  fBufferCount = 0;
  fBufferPointer[METAFILE_RECORD_LENGTH] = '\0';
  fPageNumber = 0;
  // Header :
  printFLN   ("%%!PS-Adobe-2.0");
  printFLN   ("%%%%Creator: HEPVis::SbPainterPS.");
  printFLN("%%%%CreationDate: %s",GetDate());
  printFLN("%%%%Title: %s",fFileName);
  printFLN("%%%%Pages: (atend)");
  printFLN("%%%%BoundingBox: 0 0 %d %d",
           (int)fDeviceWidth,(int)fDeviceHeight);
  printFLN("%%%%DocumentFonts: Courier-Bold");
  printFLN("%%%%DocumentPaperSizes: a4");
  printFLN("%%%%EndComments");
  // PostScript :
  putSaveStateInStream      ();
  // General :
  putInStreamF("/n {newpath} def ");
  putInStreamF("/cl {closepath} def ");
  putInStreamF("/s {stroke} def ");
  putInStreamF("/f {fill} def ");
  // Move :
  putInStreamF("/m  {moveto} def ");
  putInStreamF("/rm {rmoveto} def ");
  putInStreamF("/rl {rlineto} def ");
  // Line :
  putInStreamF("/lc {setlinecap} def ");
  putInStreamF("/lw {setlinewidth} def ");
  putInStreamF("/rgb {setrgbcolor} def ");
  putInStreamF("/ss {[] 0 setdash} def ") ;            /* style solid       */
  putInStreamF("/sd {[12 6] 0 setdash} def ");         /* style dashed      */
  putInStreamF("/so {[6 12] 0 setdash} def ");         /* style dotted      */
  putInStreamF("/sdo {[18 12 6 12] 0 setdash} def ");  /* style dash dotted */
  // Mark :
  fMarkerSize = 2.;
  putInStreamF("/ms 2. def /msi .5 def ");        /* mark size */
  putInStreamF("/cross {ms ms scale -1. -1. rm  ");
  putInStreamF("2. 2. rl 0. -2. rm -2. 2. rl msi msi scale} def ");
  putInStreamF("/plus  {ms ms scale -1. 0. rm 2. 0. rl ");
  putInStreamF("-1. 1. rm 0. -2. rl msi msi scale} def ");
  putInStreamF("/asterisk {ms ms scale -1. 0. rm 2. 0. rl -1. 1. rm ");
  putInStreamF("0. -2. rl 0. 1. rm -0.707 -0.707 rm 1.414 1.414 rl ");
  putInStreamF("0. -1.414 rm -1.414 1.414 rl msi msi scale} def ");
  putInStreamF("/triangle {ms ms scale 0. 1. rm -0.6 -1.5 rl ");
  putInStreamF("1.2 0. rl -0.6 1.5 rl msi msi scale} def ");
  // Text :
  putInStreamF("/sh {show} def ");
  putInStreamF("/df {/Courier-Bold findfont} def ");
  putInStreamF("/mf {makefont setfont} def ");
  printFLN("%%%%EndProlog");
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::closeStream(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile==NULL) return;
  putRestoreStateInStream    ();
  printFLN("%%%%Trailer");
  printFLN("%%%%Pages: %d",fPageNumber);
  printFLN("%%%%EOF");
  if(fFile!=NULL) ::fclose(fFile);
  fFile = NULL;
  STRDEL(fFileName);
  fFileName = NULL;
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putInStreamF(
 const char* aFormat 
,...
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile==NULL) return;
  va_list  args;
  va_start(args,aFormat);
  printV(aFormat,args);
  va_end(args);
  int length = (int)strlen(fBufferString);
  if(length>METAFILE_RECORD_LENGTH) {
    ::printf("SoPostScript::putInStreamF overflow\n");
    return;
  }
  int nlength = fBufferCount + length;
  if(nlength>METAFILE_RECORD_LENGTH) {
      fBufferPointer[fBufferCount] = '\0';
      if(::fprintf(fFile,"%s\n",(char*)fBufferPointer)<0) {
        ::printf("SoPostScript::putInStreamF fprintf error\n");
      }
      fBufferCount = 0;
      nlength = length;
    }
  unsigned char* pointer = fBufferPointer + fBufferCount;
  ::strcpy((char*)pointer,fBufferString);
  fBufferCount = nlength;
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::printFLN(
 const char* aFormat 
,...
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fFile==NULL) return;
  va_list args;
  va_start(args,aFormat);
  printV(aFormat,args);
  va_end(args);
/* put buffer in file */
  if(fBufferCount>0) {
    fBufferPointer[fBufferCount] = '\0';
    if(::fprintf (fFile,"%s\n",(char*)fBufferPointer)<0) {
      ::printf("SbPainterPS::printFLN fprintf error\n");
    }
    fBufferCount = 0;
  }
/* put comment in file */
  if(::fprintf (fFile,"%s\n",fBufferString)<0) {
    ::printf("SbPainterPS::printFLN fprintf error\n");
  }
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::printV(
 const char* This 
,va_list aArgs  
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#define MAX_STR    2048
  if(fBufferString==NULL) {
    fBufferString = (char*)::malloc(MAX_STR * sizeof(char));
    if(fBufferString==NULL) return;
  }
  fBufferString[MAX_STR-1]  = '\0';
  ::vsnprintf(fBufferString,MAX_STR-1, This,aArgs);
  if(fBufferString[MAX_STR-1]!='\0') {
    ::printf("SbPainterPS::printV overflow\n");
    fBufferString[0] = '\0';
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putPageScaleInStream(
 float aWidth 
,float aHeight 
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(aWidth <=0.) aWidth  = 100.;
  if(aHeight<=0.) aHeight = 100.;

  putScaleInStream (1./METAFILE_SCALE,1./METAFILE_SCALE);
  putTranslationInStream ((float)(fDeviceWidth/20.),
                          (float)(fDeviceHeight/30.));

  float scale;
  if(fDeviceWidth<=fDeviceHeight)
    scale = (aHeight<=aWidth ? 
             fDeviceWidth /aWidth  : fDeviceWidth /aHeight );
  else 
    scale = (aHeight<=aWidth ? 
             fDeviceHeight /aWidth : fDeviceHeight /aHeight );

  float xtra,ytra;
  if(fParams.portrait==1) {
    xtra = (fDeviceWidth  - scale * aWidth)/2;
    ytra = (fDeviceHeight - scale * aHeight)/2;
  } else {
    putTranslationInStream(fDeviceWidth,0.);
    putRotateInStream(90);
    xtra = (fDeviceHeight  - scale * aWidth)/2;
    ytra = (fDeviceWidth   - scale * aHeight)/2;
  }
  putTranslationInStream (xtra,ytra);

  putScaleInStream (scale,scale);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putSaveStateInStream(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("gsave ");
  fGSave++;
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putRestoreStateInStream(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("grestore ");
  fGSave--;
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putTranslationInStream(
 float aX
,float aY
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("%.2f %.2f translate ",aX,aY);
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putScaleInStream(
 float aX
,float aY
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("%.2f %.2f scale ",aX,aY);
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putBeginPageInStream(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fPageNumber++;
  printFLN("%%%%Page: %d %d",fPageNumber,fPageNumber);
  putSaveStateInStream();
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putEndPageInStream (
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("showpage ");
  putRestoreStateInStream();
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putRGB_InStream (
 float aR 
,float aG
,float aB 
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fParams.shade==Color)       
    putInStreamF("%.2f %.2f %.2f rgb ",aR,aG,aB);
  else if(fParams.shade==Grey)
    putInStreamF("%.2f setgray ",convertRGB_ToGrey(aR,aG,aB));
  else if(fParams.shade==BlackWhite)  
    putInStreamF("0. setgray ",convertRGB_ToGrey(aR,aG,aB));
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putLineWidthInStream(
 int aWidth
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fParams.lineWidth<0.) {
    if(aWidth==1) {
      putInStreamF("%.1f lw ",0.5); // For a better rendering.
    } else {
      putInStreamF("%.1f lw ",(float)(aWidth));
    }
  } else {
    putInStreamF("%.1f lw ",fParams.lineWidth); 
  }
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putMarkerSizeInStream (
 float aSize
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(aSize==fMarkerSize) return;
  fMarkerSize = aSize;
  putInStreamF("/ms %g def /msi %g def ",aSize,1./aSize);
}
/*
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putMarkerStyleInStream (
 SbMarkerStyle aStyle 
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  switch (aStyle) {
  case SbMarkerPlus:
    putInStreamF("plus ");
    break;
  case SbMarkerAsterisk:
  case SbMarkerStar:
    putInStreamF("asterisk ");
    break;
  case SbMarkerCross:
    putInStreamF("cross ");
    break;
  case SbMarkerTriangleUpLine:
    putInStreamF("triangle ");
    break;
  default:
    putLineToInStream(0.,0.);
    break;
  }
}
*/
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putBackgroundInStream (
 float aR
,float aG
,float aB
,float aWidth 
,float aHeight 
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putNewPathInStream(); 
  putMoveInStream(0.,0.);
  putLineToInStream(aWidth,0.);
  putLineToInStream(0.,aHeight);
  putLineToInStream(-aWidth,0.);
  putLineToInStream(0.,-aHeight);
  putClosePathInStream();
  if(fParams.doBack==1) {
    // Back :
    putSaveStateInStream();
    putRGB_InStream(aR,aG,aB);
    putFillInStream();       
    putRestoreStateInStream();
  }
  // Clip :
  putInStreamF("clip ");
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putFrameInStream (
 float aR
,float aG
,float aB
,float aWidth 
,float aHeight 
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putNewPathInStream(); 
  putMoveInStream(0.,0.);
  putLineToInStream(aWidth,0.);
  putLineToInStream(0.,aHeight);
  putLineToInStream(-aWidth,0.);
  putLineToInStream(0.,-aHeight);
  putClosePathInStream();
  putRGB_InStream(aR,aG,aB);
  putLineWidthInStream(1);
  putCapInStream(1);
  putInStreamF("ss ");
  putStrokeInStream();
}
//////////////////////////////////////////////////////////////////////////////
float SbPainterPS::convertRGB_ToGrey (
 float aRed
,float aGreen
,float aBlue
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return (0.3F * aRed + 0.59F * aGreen + 0.11F * aBlue);
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putRotateInStream(
 float aX                  
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("%.2f  rotate ",aX);
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putNewPathInStream(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("n ");
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putStrokeInStream(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("s ");
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putFillInStream(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("f ");
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putClosePathInStream(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("cl ");
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putCapInStream(
 int aX
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("%1d lc ",aX);
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putLineToInStream(
 float aX
,float aY
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF ("%.2f %.2f rl ",aX,aY);
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putMoveInStream(
 float aX
,float aY
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF ("%.2f %.2f m ",aX,aY);
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putCircleInStream(
 float aX                            
,float aY                            
,float aR                            
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  putInStreamF("%.2f %.2f %.2f 0 360 arc s ",aX,aY,aR);
}
/*
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putLineStyleInStream(
 SbLineStyle aStyle
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  switch(aStyle) {
  case SbLineSolid:putInStreamF("ss ") ;break;
  case SbLineDashed:putInStreamF("sd ") ;break;
  case SbLineDotted:putInStreamF("so ") ;break;
  case SbLineDashDotted:putInStreamF("sdo ");break;
  }
}
*/
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/////// Image ////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::putImageInStream (
 unsigned int aWidth
,unsigned int aHeight
,GetRGB_Function aProc
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if((aWidth<=0)||(aHeight<=0)) return;
  if(!aProc) return;

  putSaveStateInStream      ();
  putInStreamF     ("%d %d scale ", aWidth, aHeight ); 
  int status = 1;
  int      nbhex;
  unsigned int row,col,col_max;
  double   dr,dg,db;
  typedef unsigned char Uchar;
  Uchar    red,green,blue,b;
  if(fParams.shade!=0) { /*grey*/
    putInStreamF   ("/picstr %d string def ",aWidth); 
    putInStreamF   ("%d %d %d ",aWidth,aHeight,8); 
    putInStreamF   ("[ %d 0 0 -%d 0 %d ] ",aWidth,aHeight,aHeight); 
    putInStreamF   ("{ currentfile picstr readhexstring pop } " );
    printFLN ("image " );
    for ( row = 0; row < aHeight; row++ ){
      for ( col = 0; col < aWidth; col++){ 
        double    fgrey;
        Uchar     grey;
        status    = aProc(col,row,dr,dg,db)==0 ? 0 : status;
        fgrey     = ConvertRGB_ToGrey(dr,dg,db);
        grey      = (Uchar) ( 255. * fgrey);
        writeByte (grey);
      }
    }
    nbhex     = aWidth * aHeight * 2; 
    printFLN ("%%%% nbhex digit          :%d ",nbhex); 
    printFLN ("%%%% nbhex/record_length  :%d ",nbhex/METAFILE_RECORD_LENGTH); 
    printFLN ("%%%% nbhex%%record_length :%d ",nbhex%METAFILE_RECORD_LENGTH); 
  }else if(fParams.nbit==2){ 
    int       nbyte2;
    nbyte2    = (aWidth   *  3)/4;
    nbyte2   /=3;
    nbyte2   *=3;
    col_max   = (nbyte2  *  4)/3;
    /* 2 bit for r and g and b   */
    /* rgbs following each other */
    putInStreamF   ("/rgbstr %d string def ",nbyte2); 
    putInStreamF   ("%d %d %d ",col_max,aHeight,2); 
    putInStreamF   ("[ %d 0 0 -%d 0 %d ] ",col_max,aHeight,aHeight); 
    putInStreamF   ("{ currentfile rgbstr readhexstring pop } " );
    putInStreamF   ("false 3 " );
    printFLN ("colorimage " );
    for ( row = 0; row < aHeight; row++ ){
      for ( col = 0; col < col_max; col+=4){
        status  = aProc(col,row,dr,dg,db)==0 ? 0 : status;
        red     = (Uchar) ( 3. * dr);
        green   = (Uchar) ( 3. * dg);
        blue    = (Uchar) ( 3. * db);
        b       = red;
        b       = (b<<2)+green;
        b       = (b<<2)+blue;
        status  = aProc(col+1,row,dr,dg,db)==0 ? 0 : status;
        red     = (Uchar) ( 3. * dr);
        green   = (Uchar) ( 3. * dg);
        blue    = (Uchar) ( 3. * db);
        b     = (b<<2)+red;
        writeByte (b);
        
        b       = green;
        b       = (b<<2)+blue;
        status  = aProc(col+2,row,dr,dg,db)==0 ? 0 : status;
        red     = (Uchar) ( 3. * dr);
        green   = (Uchar) ( 3. * dg);
        blue    = (Uchar) ( 3. * db);
        b     = (b<<2)+red;
        b     = (b<<2)+green;
        writeByte (b);
        
        b       = blue;
        status  = aProc(col+3,row,dr,dg,db)==0 ? 0 : status;
        red     = (Uchar) ( 3. * dr);
        green   = (Uchar) ( 3. * dg);
        blue    = (Uchar) ( 3. * db);
        b     = (b<<2)+red;
        b     = (b<<2)+green;
        b     = (b<<2)+blue;
        writeByte (b);
      }
    }
  }else if(fParams.nbit==4){ 
    int       nbyte4;
    nbyte4    = (aWidth  * 3)/2;
    nbyte4   /=3;
    nbyte4   *=3;
    col_max   = (nbyte4 * 2)/3;
    /* 4 bit for r and g and b   */
    /* rgbs following each other */
    putInStreamF   ("/rgbstr %d string def ",nbyte4); 
    putInStreamF   ("%d %d %d ",col_max,aHeight,4); 
    putInStreamF   ("[ %d 0 0 -%d 0 %d ] ",col_max,aHeight,aHeight); 
    putInStreamF   ("{ currentfile rgbstr readhexstring pop } " );
    putInStreamF   ("false 3 " );
    printFLN ("colorimage " );
    for ( row = 0; row < aHeight; row++ ){
      for ( col = 0; col < col_max; col+=2){
        status  = aProc(col,row,dr,dg,db)==0 ? 0 : status;
        red     = (Uchar) ( 15. * dr);
        green   = (Uchar) ( 15. * dg);
        putInStreamF ("%x%x",red,green);
        blue    = (Uchar) ( 15. * db);
        
        status  = aProc(col+1,row,dr,dg,db)==0 ? 0 : status;
        red     = (Uchar) ( 15. * dr);
        putInStreamF ("%x%x",blue,red);
        green   = (Uchar) ( 15. * dg);
        blue    = (Uchar) ( 15. * db);
        putInStreamF ("%x%x",green,blue);
      }
    }
  }else{ 
    int       nbyte8;
    nbyte8    = aWidth   * 3;
    /* 8 bit for r and g and b   */
    putInStreamF   ("/rgbstr %d string def ",nbyte8); 
    putInStreamF   ("%d %d %d ",aWidth,aHeight,8); 
    putInStreamF   ("[ %d 0 0 -%d 0 %d ] ",aWidth,aHeight,aHeight); 
    putInStreamF   ("{ currentfile rgbstr readhexstring pop } " );
    putInStreamF   ("false 3 " );
    printFLN   ("colorimage " );
    for ( row = 0; row < aHeight; row++ ){
      for ( col = 0; col < aWidth; col++){
        status     = aProc(col,row,dr,dg,db)==0 ? 0 : status;
        red        = (Uchar) ( 255. * dr);
        writeByte (red);
        green      = (Uchar) ( 255. * dg);
        writeByte (green);
        blue       = (Uchar) ( 255. * db);
        writeByte (blue);
      }
    }
  }
  if(status==0) 
    ::printf("SbPainterPS::putImageInStream: problem to retrieve some pixel rgb.\n");
  putRestoreStateInStream();
}
//////////////////////////////////////////////////////////////////////////////
void SbPainterPS::writeByte (
 unsigned char a_byte
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  unsigned char h = a_byte / 16;
  unsigned char l = a_byte % 16;
  putInStreamF ("%x%x",h,l);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
char* GetDate (
)
//////////////////////////////////////////////////////////////////////////////
// Return local date.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  time_t d;
  time(&d);
  char* string = ctime(&d);
  string[24] = '\0';
  return string;
}
//////////////////////////////////////////////////////////////////////////////
double ConvertRGB_ToGrey(
 double a_red
,double a_green
,double a_blue
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return (0.30 * a_red + 0.59 * a_green + 0.11 * a_blue);
}
