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
// $Id: SoTubs.h,v 1.2 2004/06/14 09:27:38 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoTubs                                                  */
/* Description:      Represents the G4Tubs Geant Geometry entity             */
/* Author:           Joe Boudreau Nov 11 1996                                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef HEPVis_SoTubs_h
#define HEPVis_SoTubs_h

#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFNode.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/nodes/SoShape.h>

class SoSFNode;

//! SoTubs - Inventor version of the G4Tubs Geant Geometry entity
/*!
 * Node:	     SoTubs
 * 
 * Description:      The Inventor version of the G4Tubs Geant Geometry entity
 * 
 * Author:	     Joe Boudreau Nov 11 1996
 * 
 * class G4Tubs
 * 
 * A tube or tube segment with curved sides parallel to
 * the z-axis. The tube has a specified half-length along
 * the z axis, about which it is centred, and a given
 * minimum and maximum radius. A minimum radius of 0
 * signifies a filled tube /cylinder. The tube segment is
 * specified by starting and delta
 * angles for phi, with 0 being the +x axis, PI/2
 * the +y axis. A delta angle of 2PI signifies a
 * complete, unsegmented tube/cylinder
 *
 * Always use Inventor Fields. This allows Inventor to detect a change to
 * the data field and take the appropriate action; e.g., redraw the scene.
*/

#define SoTubs Geant4_SoTubs

class SoTubs:public SoShape {

  // The following is required:
  SO_NODE_HEADER(SoTubs);
 
public:

  //
  //! Inside radius of the tube
  //
  SoSFFloat pRMin;   
  //
  //! Outside radius of the tube
  //
  SoSFFloat pRMax;
  //
  //! Half-length in Z
  //
  SoSFFloat pDz;
  //
  //! Starting angle, in radians
  //
  SoSFFloat pSPhi;
  //
  //! Delta-angle, in radians
  //
  SoSFFloat pDPhi;             
  //
  //! Alternate rep - required
  //
  SoSFNode  alternateRep;

  //
  //! Constructor, required
  //
  SoTubs();

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
  virtual ~SoTubs();

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
