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
#ifndef HEPVis_SoMarkerSet_h
#define HEPVis_SoMarkerSet_h 

#include <Inventor/nodes/SoPointSet.h>

#include <Inventor/fields/SoMFInt32.h>

#define HEPVis_SoMarkerSet Geant4_HEPVis_SoMarkerSet

class HEPVis_SoMarkerSet : public SoPointSet {
  SO_NODE_HEADER(HEPVis_SoMarkerSet);
public:
  SoMFInt32 markerIndex ;
  HEPVis_SoMarkerSet();

  enum MarkerType {
    //************ 5x5 markers
    PLUS_5_5 = 0,
    ASTERISK_5_5,
    CROSS_5_5,
    STAR_5_5,
    CIRCLE_LINE_5_5,
    CIRCLE_FILLED_5_5,
    TRIANGLE_UP_LINE_5_5,
    TRIANGLE_UP_FILLED_5_5,
    TRIANGLE_DOWN_LINE_5_5,
    TRIANGLE_DOWN_FILLED_5_5,
    DAVID_STAR_LINE_5_5 = 10,
    DAVID_STAR_FILLED_5_5,
    SWISS_CROSS_LINE_5_5,
    SWISS_CROSS_FILLED_5_5,
    DIAMOND_LINE_5_5,
    DIAMOND_FILLED_5_5,
    SQUARE_LINE_5_5,
    SQUARE_FILLED_5_5 = 17,

    //************ 7x7 markers
    PLUS_7_7 = 18,
    ASTERISK_7_7,
    CROSS_7_7,
    STAR_7_7,
    CIRCLE_LINE_7_7,
    CIRCLE_FILLED_7_7,
    TRIANGLE_UP_LINE_7_7,
    TRIANGLE_UP_FILLED_7_7,
    TRIANGLE_DOWN_LINE_7_7,
    TRIANGLE_DOWN_FILLED_7_7,
    DAVID_STAR_LINE_7_7,
    DAVID_STAR_FILLED_7_7,
    SWISS_CROSS_LINE_7_7 = 30,
    SWISS_CROSS_FILLED_7_7,
    DIAMOND_LINE_7_7,
    DIAMOND_FILLED_7_7,
    SQUARE_LINE_7_7,
    SQUARE_FILLED_7_7 = 35,

    //************ 9x9 markers
    PLUS_9_9 = 36,
    ASTERISK_9_9,
    CROSS_9_9,
    STAR_9_9,
    CIRCLE_LINE_9_9,
    CIRCLE_FILLED_9_9,
    TRIANGLE_UP_LINE_9_9,
    TRIANGLE_UP_FILLED_9_9,
    TRIANGLE_DOWN_LINE_9_9,
    TRIANGLE_DOWN_FILLED_9_9,
    DAVID_STAR_LINE_9_9,
    DAVID_STAR_FILLED_9_9,
    SWISS_CROSS_LINE_9_9 = 30,
    SWISS_CROSS_FILLED_9_9,
    DIAMOND_LINE_9_9,
    DIAMOND_FILLED_9_9,
    SQUARE_LINE_9_9,
    SQUARE_FILLED_9_9 = 53
  } ;

public: /*SoEXTENDER*/
  virtual void GLRender(SoGLRenderAction*);
public: /*SoINTERNAL*/
  static void initClass();
protected:
  virtual ~HEPVis_SoMarkerSet();
};

#endif
