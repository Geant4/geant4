/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoBox                                                 */
/* Description:      Represents the G4Box Geant Geometry entity              */
/* Author:           Joe Boudreau Nov 11 1996                                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef HEPVis_SoBox_h
#define HEPVis_SoBox_h

#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFNode.h>
#include <Inventor/nodes/SoShape.h>

class SoSFNode;

//! SoBox - Inventor version of the G4Box Geant Geometry entity
/*!
 * Node:	     SoBox
 * 
 * Description:      The Inventor version of the G4Box Geant Geometry entity
 * 
 * Author:	     Joe Boudreau Nov 11 1996
 * 
 * The documentation from Geant says:
 *
 * class G4Box
 *
 * A Box is a cuboid of given half lengths dx,dy,dz. The Box is
 * centred on the origin with sides parallel to the x/y/z axes.
 *
 * Always use Inventor Fields. This allows Inventor to detect a change to
 * the data field and take the appropriate action; e.g., redraw the scene.
*/

#define SoBox Geant4_SoBox

class SoBox:public SoShape {

  // The following is required:
  SO_NODE_HEADER(SoBox);
 
public:

  //
  //! Half-length along X
  //
  SoSFFloat fDx;   
  //
  //! Half-length along Y
  //
  SoSFFloat fDy;
  //
  //! Half-length along Z
  //
  SoSFFloat fDz;
  //
  //! Alternate rep - for use by users without HEPVis shared objects
  //
  SoSFNode  alternateRep;

  //
  //! Constructor, required
  //
  SoBox();

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
  virtual ~SoBox();

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

};

#endif
