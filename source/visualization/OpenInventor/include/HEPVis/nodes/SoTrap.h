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
//
// $Id: SoTrap.h 66373 2012-12-18 09:41:34Z gcosmo $
//
/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoTrap                                                  */
/* Description:      Represents the G4Trap Geant Geometry entity             */
/* Author:           Joe Boudreau Nov 11 1996                                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef HEPVis_SoTrap_h
#define HEPVis_SoTrap_h

// Inheritance :
#include <Inventor/nodes/SoShape.h>

#include <Inventor/fields/SoSFNode.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFFloat.h>

//adding fields 
#ifdef WIN32
#include <SoWinEnterScope.h>
#endif


// Defining Class
#ifdef WIN32
#include <SoWinLeaveScope.h>
#endif

class SoSFNode;

//! SoTrap - Inventor version of the G4Trap Geant Geometry entity
/*! Node:          SoTrap
 *
 * Description:	   Inventor version of the G4Trap Geant Geometry entity
 *
 * Author:	   Joe Boudreau Nov 11 1996
 *
 * A G4Trap is a general trapezoid: The faces perpendicular to the z planes
 * are tapezia, and their centres are not necessarily on a line parallel to
 * the z axis.
 *
 * Note that of the 11 parameters desribed below, only 9 are really
 * independent - a check for planarity is made in the calculation of the
 * equation for each plane. If the planes are not parallel, a call to
 * G4Exception is made.
 *
 * Always use Inventor Fields. This allows Inventor to detect a change to
 * the data field and take the appropriate action; e.g., redraw the scene.
 *
*/

#define SoTrap Geant4_SoTrap

class SoTrap:public SoShape {

  // The following is required:
  SO_NODE_HEADER(SoTrap);
 
public:
#if defined(WIN32) && defined(BUILDING_DLL)
   // When building the node as a DLL under Win32 we must explicitly
   // declare this entry point as visible outside the DLL.  The macro
   // BUILDING_DLL is defined in the node's source file.
      _declspec(dllexport)
#endif

  //
  //! half-length along Z
  //
  SoSFFloat pDz;                  
  //
  //! Polar angle of the line joining the centres of the faces at -/+pDz
  //
  SoSFFloat pTheta;
  //
  //! Azimuthal angle of the line joing the centre of the face at -pDz 
  //! to the centre of the face at +pDz
  //
  SoSFFloat pPhi;
  //
  //! Half-length along y of the face at -pDz
  //
  SoSFFloat pDy1;
  //
  //! Half-length along x of the side at y=-pDy1 of the face at -pDz
  //
  SoSFFloat pDx1;
  //
  //! Half-length along x of the side at y=+pDy1 of the face at -pDz
  //
  SoSFFloat pDx2;
  //
  //! Half-length along y of the face at +pDz
  //
  SoSFFloat pDy2;
  //
  //! Half-length along x of the side at y=-pDy2 of the face at +pDz
  //
  SoSFFloat pDx3;
  //
  //! Half-length along x of the side at y=+pDy2 of the face at +pDz
  //
  SoSFFloat pDx4;
  //
  //! Angle with respect to the y axis from the centre of the side at 
  //! y=-pDy1 to the centre at y=+pDy1 of the face at -pDz
  //
  SoSFFloat pAlp1;
  //
  //! Angle with respect to the y axis from the centre of the side at 
  //! y=-pDy2 to the centre at y=+pDy2 of the face at +pDz
  //
  SoSFFloat pAlp2;
 
  //
  //! Alternate rep - required
  //
  SoSFNode  alternateRep;

  //
  //! Constructor, required
  //
  SoTrap();

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

  //
  //! Destructor, required
  //
  virtual ~SoTrap();

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

#ifdef WIN32
#include <SoWinEnterScope.h>
#endif

#endif
