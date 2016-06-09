//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// History :
// - 09/12/98, G.Barrand : creation.

#ifndef HEPVis_SbPainter_h
#define HEPVis_SbPainter_h 

#include <Inventor/SbBasic.h>

/*
#include <HEPVis/SbStyle.h>

typedef enum {
  SbPrimitivePoints,
  SbPrimitiveLineStrip,
  SbPrimitiveLineLoop,
  SbPrimitiveLines,
  SbPrimitivePolygon
} SbPrimitiveType;

class SbPainterContext {
public:
  float fRed,fGreen,fBlue;
  int fLineWidth;
  SbLineStyle fLineStyle;
  int fMarkerSize;
  SbMarkerStyle fMarkerStyle;
  SbAreaStyle fAreaStyle;
  int fPolygonMode;
};
*/

#define SbPainter Geant4_SbPainter

class SbPainter {
public:
  SbPainter();
  virtual ~SbPainter();
  void setWindowSize(int,int);
  void enableEdges(SbBool);
  virtual void beginTraversal() = 0;
  virtual void clearColorBuffer(float,float,float) = 0;
  //virtual void drawPrimitive(SbPrimitiveType,
  //                          int,float*,float*,float*,
   //                        const SbPainterContext&) = 0;
  virtual void endTraversal() = 0;
protected:
  int fWindowWidth,fWindowHeight;
  float fRed,fGreen,fBlue;
  SbBool fEdges;
};

#endif
