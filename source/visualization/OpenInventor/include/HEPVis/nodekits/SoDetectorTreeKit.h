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
// $Id: SoDetectorTreeKit.h,v 1.3 2004/11/22 22:57:01 gbarrand Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoDetectorTreeKit                                       */
/* Description:      Easy way of browsing through a tree of detectors        */
/* Author:           Joe Boudreau Nov 11 1996                                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef HEPVis_SoDetectorTreeKit_h
#define HEPVis_SoDetectorTreeKit_h

#include <Inventor/nodekits/SoBaseKit.h>
#include <Inventor/nodes/SoEventCallback.h>

class SoSFNode;
class SoHandleEventAction;
class SoSeparator;

#define SoDetectorTreeKit Geant4_SoDetectorTreeKit


class SoDetectorTreeKit:public SoBaseKit {

  // The following is required:
  SO_KIT_HEADER(SoDetectorTreeKit);
  ////////////////////////////////////////////
  SoSFNode alternateRep;
  ////////////////////////////////////////////
  SO_KIT_CATALOG_ENTRY_HEADER(callbackList);
  SO_KIT_CATALOG_ENTRY_HEADER(topSeparator);
    SO_KIT_CATALOG_ENTRY_HEADER(pickStyle);
    SO_KIT_CATALOG_ENTRY_HEADER(appearance);
    SO_KIT_CATALOG_ENTRY_HEADER(units);
    SO_KIT_CATALOG_ENTRY_HEADER(transform);
    SO_KIT_CATALOG_ENTRY_HEADER(texture2Transform);
    SO_KIT_CATALOG_ENTRY_HEADER(childList);
      SO_KIT_CATALOG_ENTRY_HEADER(previewSeparator);
      SO_KIT_CATALOG_ENTRY_HEADER(fullSeparator);
   

public:

  // Constructor, required
  SoDetectorTreeKit();

  // This is required
  virtual SbBool affectsState() const;

  // Class Initializer, required
  static void initClass();

  // Turn the preview on or off
  virtual void setPreview(SbBool Flag);

  // Return the preview state
  virtual SbBool getPreview() const;

  // The SoSwitch::whichChild = SO_SWITCH_ALL
  void setPreviewAndFull();

  // Return the preview Separator
  virtual SoSeparator *getPreviewSeparator() const;

  // Return the full Separator
  virtual SoSeparator *getFullSeparator() const;

  // Generate AlternateRep, required.  Generating an alternate representation
  // must be done upon users request.  It allows an Inventor program to read
  // back the file without requiring *this* code to be dynamically linked. 
  // If the users expects that *this* code will be dynamically linked, he
  // need not invoke this method.  
  virtual void generateAlternateRep();

  // We better be able to clear it, too!
  virtual void clearAlternateRep();

protected:

  // Destructor.
  virtual ~SoDetectorTreeKit();

private: 

  // This is needed as well
  void createInitialTree();

  // This is the callback function that will be 
  // added to the callback list
  static void   expand   (void *userData, SoEventCallback *eventCB);
  static void   contract (void *userData, SoEventCallback *eventCB);
};

#endif
