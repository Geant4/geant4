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
//
// $Id: SoTrd.h,v 1.2 2004/06/14 09:27:37 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoTrd                                                   */
/* Description:      Represents the G4Trd Geant Geometry entity              */
/* Author:           Joe Boudreau Nov 11 1996                                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef HEPVis_SoTrd_h
#define HEPVis_SoTrd_h

#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFNode.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/nodes/SoShape.h>

class SoSFNode;

//! SoTrd - Inventor version of the G4Trd Geant Geometry entity
/*!
 * Node:            SoTrd
 * 
 * Description:     Inventor version of the G4Trd Geant Geometry entity
 * 
 * Author:	    Joe Boudreau Nov 11 1996
 * 
 * A Trd is a trapezoid with the x and y dimensions varying along z
 *
 * Always use Inventor Fields. This allows Inventor to detect a change to
 * the data field and take the appropriate action; e.g., redraw the scene.
 *
 */

#define SoTrd Geant4_SoTrd

class SoTrd:public SoShape {

  //
  //! This is required
  //
  SO_NODE_HEADER(SoTrd);
 
public:

  //
  //! half-length of x, at -fDz
  //
  SoSFFloat fDx1;                  

  //
  //! half-length of x, at +fDz
  //
  SoSFFloat fDx2;

  //
  //! half-length of y, at -fDz
  //
  SoSFFloat fDy1;

  //
  //! half-length of y, at +fDz
  //
  SoSFFloat fDy2;

  //
  //! half-length along Z
  //
  SoSFFloat fDz;

  //
  //! Alternate rep - required
  //
  SoSFNode  alternateRep;

  //
  //! Constructor, required
  //
  SoTrd();

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
  virtual ~SoTrd();

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
