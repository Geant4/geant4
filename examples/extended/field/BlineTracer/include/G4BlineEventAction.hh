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
/// \file field/BlineTracer/include/G4BlineEventAction.hh
/// \brief Definition of the G4BlineEventAction class
//
//
//
// 
// --------------------------------------------------------------------
//
// G4BlineEventAction
//
// Class description:
//
// Defines the EventAction used during tracing of magnetic field
// lines. During the EndOfEventAction() it stores magneticfield lines
// as polyline and polymarker objects.
// These polyline and polymarker objects can be drawn by using the
// DrawFieldLines() method.

// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------
#ifndef G4BlineEventAction_h
#define G4BlineEventAction_h 1 

#include "G4UserEventAction.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

class G4Event;
class G4Polyline;
class G4Polymarker;
class G4BlineTracer;

class G4BlineEventAction : public G4UserEventAction
{
  public:  // with description

    G4BlineEventAction(G4BlineTracer* aBlineTool);
    virtual ~G4BlineEventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void DrawFieldLines(G4double zoom, G4double theta, G4double phi); 
    void ResetVectorObjectToBeDrawn();

  public:  // with description

    inline void SetDrawColour(G4Colour aColour) { fDrawColour = aColour; }
    inline void SetDrawBline(G4bool aBool) { fDrawBline=aBool; }
    inline void SetDrawPoints(G4bool aBool) { fDrawPoints=aBool; }
    inline void SetPointSize(G4double aVal) { fPointSize=aVal; }
    inline G4bool GetDrawBline() { return fDrawBline; }

  public:  // without description

    // inline void SetDrawLineWidth(G4double aVal)
    //   { fTrajectoryVisAttributes.SetLineWidth(aVal); }
    // inline void SetDrawLineStyle(G4VisAttributes::LineStyle aStyle)
    //   { fTrajectoryVisAttributes.SetLineStyle(aStyle); }
      // Future implementation...

  private:
  
   G4BlineTracer* fBlineTool;   
   G4Colour fDrawColour;
   G4bool fDrawBline;
   G4bool fDrawPoints;
   G4double fPointSize;
   std::vector<G4VisAttributes*> fTrajectoryVisAttributes;
   std::vector<G4Polyline> fTrajectoryPolyline;
   std::vector<G4Polymarker> fTrajectoryPoints;
};

#endif
