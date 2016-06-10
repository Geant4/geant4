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
// $Id: G4GeometryCellStep.hh 66356 2012-12-18 09:02:32Z gcosmo $
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
