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
// $Id: G4FaceOuterBoundCreator.hh,v 1.5 2002-11-21 16:49:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4FaceOuterBoundCreator
//
// Class description:
//
//

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------
#ifndef G4FACEOUTERBOUNDCREATOR_HH
#define G4FACEOUTERBOUNDCREATOR_HH

#include "G4FaceBoundCreator.hh"

class G4FaceOuterBoundCreator: public G4FaceBoundCreator
{
  public:

  // Constructor & destructor

    G4FaceOuterBoundCreator();
    ~G4FaceOuterBoundCreator();

  // Member functions

  void CreateSTEPGeometry(void* G4obj);
  const char* Name() const { return "Face_Outer_Bound"; }
  static G4FaceOuterBoundCreator GetInstance();

  // Members

  private:

    static G4FaceOuterBoundCreator csc;
};

#endif
