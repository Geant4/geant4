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
// $Id: G4VisTrajContext.hh,v 1.1 2006-05-02 20:47:40 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay May 2006
//
#ifndef G4VISTRAJCONTEXT_HH
#define G4VISTRAJCONTEXT_HH

#include "G4Colour.hh"
#include "G4Polymarker.hh"

class G4VisTrajContext {

public:

  // Default configuration
  G4VisTrajContext(const G4String& name = "Unspecified");

  // Destructor
  virtual ~G4VisTrajContext();

  G4String Name() const;
  
  void SetVisible(const G4bool& visible); 

  // Print configuration
  void Print(std::ostream& ostr) const;

  // Line
  void SetLineColour(const G4Colour& colour);
  G4Colour GetLineColour() const;

  void SetDrawLine(const G4bool& draw); 
  G4bool GetDrawLine() const;

  void SetLineVisible(const G4bool& visible);
  G4bool GetLineVisible() const;

  // Auxiliary points
  void SetDrawAuxPts(const G4bool& draw);
  G4bool GetDrawAuxPts() const;

  void SetAuxPtsType(const G4Polymarker::MarkerType& marker);
  G4Polymarker::MarkerType GetAuxPtsType() const;

  void SetAuxPtsSize(const G4double& size);
  G4double GetAuxPtsSize() const;

  void SetAuxPtsFillStyle(const G4VMarker::FillStyle& style);
  G4VMarker::FillStyle GetAuxPtsFillStyle() const;

  void SetAuxPtsColour(const G4Colour& colour);
  G4Colour GetAuxPtsColour() const;

  void SetAuxPtsVisible(const G4bool& visible);
  G4bool GetAuxPtsVisible() const;

  // Step points
  void SetDrawStepPts(const G4bool& draw);
  G4bool GetDrawStepPts() const;

  void SetStepPtsType(const G4Polymarker::MarkerType& marker);
  G4Polymarker::MarkerType GetStepPtsType() const;

  void SetStepPtsSize(const G4double& size);
  G4double GetStepPtsSize() const;

  void SetStepPtsFillStyle(const G4VMarker::FillStyle& style);
  G4VMarker::FillStyle GetStepPtsFillStyle() const;

  void SetStepPtsColour(const G4Colour& colour);
  G4Colour GetStepPtsColour() const;

  void SetStepPtsVisible(const G4bool& visible);
  G4bool GetStepPtsVisible() const;

private:
  
  // Data members
  G4String fName;

  // Line data
  G4Colour fLineColour;
  G4bool fLineVisible;
  G4bool fDrawLine;

  // Auxiliary point data
  G4bool fDrawAuxPts;
  G4Polymarker::MarkerType fAuxPtsType;
  G4double fAuxPtsSize;
  G4VMarker::FillStyle fAuxPtsFillStyle;
  G4Colour fAuxPtsColour;
  G4bool fAuxPtsVisible;

  // Step point data
  G4bool fDrawStepPts;
  G4Polymarker::MarkerType fStepPtsType;
  G4double fStepPtsSize;
  G4VMarker::FillStyle fStepPtsFillStyle;
  G4Colour fStepPtsColour;
  G4bool fStepPtsVisible;

};

#include "G4VisTrajContext.icc"

#endif

