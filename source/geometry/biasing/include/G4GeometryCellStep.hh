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
// $Id: G4GeometryCellStep.hh,v 1.1 2002-10-22 13:18:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4GeometryCellStep
//
// Class description:
//
// This class serves to address the "cell" a track previously 
// touched and a "cell" a track is currently in. It is used 
// for scoring and importance sampling in the "mass" geometry as well 
// as in a "parallel" geometry. 
// The "cell" information is available with the GetPreGeometryCell and 
// the GetPostGeometryCell functions.
// The GetCrossBoundary finction returns true in case the step
// crosses a boundary in the geometry this G4GeometryCellStep 
// refers to.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4GeometryCellStep_hh
#define G4GeometryCellStep_hh G4GeometryCellStep_hh

#include "G4GeometryCell.hh"

class G4GeometryCellStep
{
public: // with description
  G4GeometryCellStep(const G4GeometryCell &preCell, 
	  const G4GeometryCell &postCell);
    // initialise pre and post G4GeometryCell 

  ~G4GeometryCellStep();

  const G4GeometryCell &GetPreGeometryCell() const; 
    // addressing the  "cell" the track previously touched 

  const G4GeometryCell &GetPostGeometryCell() const; 
    // addressing the current "cell"  

  G4bool GetCrossBoundary() const; 
    // returns true if step crosses boundary of the geometry it refers to

  // the following functions are used by the scoring and importance
  // system to set the cell information.
  void SetPreGeometryCell(const G4GeometryCell &preCell); 
    
  void SetPostGeometryCell(const G4GeometryCell &postCell);
  
  void SetCrossBoundary(G4bool b);

private:
  G4GeometryCell fPreGeometryCell;
  G4GeometryCell fPostGeometryCell;  
  G4bool fCrossBoundary;

};



inline void G4GeometryCellStep::SetPreGeometryCell(const G4GeometryCell &preCell) 
{
  fPreGeometryCell = preCell;
}

inline void G4GeometryCellStep::SetPostGeometryCell(const G4GeometryCell &postCell)  
{
  fPostGeometryCell = postCell;
}
  
inline void G4GeometryCellStep::SetCrossBoundary(G4bool b) 
{
  fCrossBoundary = b;
}

inline const G4GeometryCell &G4GeometryCellStep::GetPreGeometryCell() const 
{
  return fPreGeometryCell;
}

inline const G4GeometryCell &G4GeometryCellStep::GetPostGeometryCell() const 
{
  return fPostGeometryCell;
}
  
inline G4bool G4GeometryCellStep::GetCrossBoundary() const 
{
  return fCrossBoundary;
}


#endif
