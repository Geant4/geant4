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
