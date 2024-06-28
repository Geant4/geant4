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
