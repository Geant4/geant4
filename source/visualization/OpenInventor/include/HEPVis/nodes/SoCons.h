/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoCons                                                */
/* Description:      Represents the G4Cons Geant Geometry entity             */
/* Author:           Joe Boudreau Nov 11 1996                                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef HEPVis_SoCons_h
#define HEPVis_SoCons_h

#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFNode.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/nodes/SoShape.h>

class SoSFNode;
//! SoCons - Inventor version of the G4Cons Geant Geometry entity
/*!
 * Node:	     SoCons
 *
 * Description:      The Inventor version of the G4Cons Geant Geometry entity
 *
 * Author:	     Joe Boudreau Nov 11 1996
 *
 * The documentation from Geant says:
 *
 * A G4Cons is, in the general case, a Phi segment of a cone, with half-length
 * fDz, inner and outer radii specified at -fDz and +fDz. The Phi segment is
 * described by a starting fSPhi angle, and the +fDPhi delta angle for the shape.
 * If the delta angle is >=2*M_PI, the shape is treated as continuous in Phi
 *
 * Member Data:
 *
 *      fRmin1  inside radius at  -fDz
 *      fRmin2  inside radius at  +fDz
 *      fRmax1  outside radius at -fDz
 *      fRmax2  outside radius at +fDz
 *      fDz     half length in z
 *
 *      fSPhi   starting angle of the segment in radians
 *      fDPhi   delta angle of the segment in radians
*/

#define SoCons Geant4_SoCons

class SoCons:public SoShape {

  // The following is required:
  SO_NODE_HEADER(SoCons);
 
public:

  //
  //! Inside radius at -fDz
  //
  SoSFFloat fRmin1;   
  //
  //! Inside radius at +fDz
  //
  SoSFFloat fRmin2;   
  //
  //! Outside radius at -fDz
  //
  SoSFFloat fRmax1;
  //
  //! Outside radius at +fDz
  //
  SoSFFloat fRmax2;
  //
  //! Half-length along Z
  //
  SoSFFloat fDz;
  //
  //! Starting angle, in radians
  //
  SoSFFloat fSPhi;
  //
  //! Delta-angle, in radians
  //
  SoSFFloat fDPhi;             
  //
  //! An Inventor option - slightly better render, worse performance
  //
  SoSFBool  smoothDraw;
  //
  //! Alternate rep required - for use by users without HEPVis shared objects
  //
  SoSFNode  alternateRep;

  //
  //! Constructor, required
  //
  SoCons();

  //
  //! Class Initializer, required
  //
  static void initClass();

  //
  //! Generate AlternateRep, required.  Generating an alternate representation
  //! must be done upon users request.  It allows an Inventor program to read
  //! back the file without requiring *this* code to be dynamically linked. 
  //! If the users expects that *this* code will be dynamically linked, he
  //! need not invoke this method.  
  //
  virtual void generateAlternateRep();

  //
  //! We better be able to clear it, too!
  //
  virtual void clearAlternateRep();

protected:

  //
  //! compute bounding Box, required
  //
  virtual void computeBBox(SoAction *action, SbBox3f &box, SbVec3f &center );

  //
  //! Generate Primitives, required
  //
  virtual void generatePrimitives(SoAction *action);

  //
  //! GetChildList, required whenever the class has hidden children
  //
  virtual SoChildList *getChildren() const;


protected:
  //
  //! Destructor, required
  //
  virtual ~SoCons();

private: 

  //
  //! Generate Children. Used to create the hidden children. Required whenever
  //! the node has hidden children.  
  //
  void generateChildren();  

  //
  //! Used to modify hidden children when a data field is changed. Required 
  //! whenever the class has hidden children. 
  //
  void updateChildren();

  //
  //! ChildList. Required whenever the class has hidden children.  
  //
  SoChildList *children;

  //
  //! help with trigonometry.  increments sines an cosines by an angle.
  //
  void inc(double & sinPhi, double & cosPhi, double sinDeltaPhi, double cosDeltaPhi) const {
    double oldSin=sinPhi,oldCos=cosPhi;
    sinPhi = oldSin*cosDeltaPhi+oldCos*sinDeltaPhi;
    cosPhi = oldCos*cosDeltaPhi-oldSin*sinDeltaPhi;    
  }
};

#endif
