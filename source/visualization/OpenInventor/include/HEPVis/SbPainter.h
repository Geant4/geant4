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
