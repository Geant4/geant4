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
// History :
// - 25/09/98, G.Barrand : creation.

#ifndef HEPVis_SbPainterPS_h
#define HEPVis_SbPainterPS_h 

// Inheritance :
#include <HEPVis/SbPainter.h>

#include <stdarg.h>

#include <Inventor/SbViewportRegion.h>

#define SbPainterPS Geant4_SbPainterPS

class SbPainterPS : public SbPainter {
public:
  SbPainterPS();
  ~SbPainterPS();
  void beginTraversal();
  void clearColorBuffer(float,float,float);
  //void drawPrimitive(SbPrimitiveType,
  //                   int,float*,float*,float*,
  //                   const SbPainterContext&);
  void endTraversal();
public:
  void setColorScheme(int);
  void setOrientation(int);
  void setBackgroundDrawn(int);
  void setBitsPerPixel(int);
  void setLineWidth(int);
  void* getStream();
  void setFileName(const char*);
  const char* getFileName() const;
  void openFileForWriting(const char*);
  void closeStream();
  void putPageScaleInStream(float,float);
  void putSaveStateInStream();
  void putRestoreStateInStream();
  void putTranslationInStream(float,float);
  void putScaleInStream(float,float);
  void putBeginPageInStream();
  void putEndPageInStream();
  void putRGB_InStream(float,float,float);
  void putMarkerSizeInStream(float);
  //void putMarkerStyleInStream(SbMarkerStyle);
  void putBackgroundInStream(float,float,float,float,float);
  void putFrameInStream(float,float,float,float,float);
  void putRotateInStream(float);
  void putNewPathInStream();
  void putStrokeInStream();
  void putFillInStream();
  void putClosePathInStream();
  void putCapInStream(int);
  void putLineToInStream(float,float);
  void putMoveInStream(float,float);
  void putCircleInStream(float,float,float);
  void putLineWidthInStream(int);
  //void putLineStyleInStream(SbLineStyle);
  typedef int(*GetRGB_Function)(unsigned int,unsigned int,
                                double&,double&,double&);
  void putImageInStream(unsigned int,unsigned int,GetRGB_Function);
private:
  //void drawPolygon(int,float*,float*,float,float,float,const SbAreaStyle&);
  //void drawLines(int,float*,float*,float,float,float,const SbLineStyle&,int);
  //void drawMarkers(int,float*,float*,
  //                 float,float,float,const SbMarkerStyle&,int);
  enum ColorScheme {
    Color = 0,
    Grey = 1,
    BlackWhite = 2
  };
  struct {
    int shade;
    int portrait;
    int nbit;
    int doBack;
    float lineWidth;
  } fParams;
  float fDeviceWidth;
  float fDeviceHeight;
  int fPageNumber;
  float fMarkerSize;
  FILE* fFile;
  char* fFileName;
  int fGSave;
  int fBufferCount;
  char* fBufferString;
#define METAFILE_RECORD_LENGTH  80
  unsigned char fBufferPointer[METAFILE_RECORD_LENGTH+1];
private:
  void putInStreamF(const char*,...);
  void printFLN(const char*,...);
  void printV(const char*,va_list);
  float convertRGB_ToGrey(float,float,float);
  void writeByte(unsigned char);
};

#endif





