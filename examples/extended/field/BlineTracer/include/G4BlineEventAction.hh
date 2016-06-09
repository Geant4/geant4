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
// $Id: G4BlineEventAction.hh,v 1.1 2003/11/25 09:29:46 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
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
    ~G4BlineEventAction();

    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

    void DrawFieldLines(G4double zoom, G4double theta, G4double phi); 
    void ResetVectorObjectToBeDrawn();

  public:  // with description

    inline void SetDrawColour(G4Colour aColour) { DrawColour = aColour; }
    inline void SetDrawBline(G4bool aBool) { DrawBline=aBool; }
    inline void SetDrawPoints(G4bool aBool) { DrawPoints=aBool; }
    inline void SetPointSize(G4double aVal) { PointSize=aVal; }
    inline G4bool GetDrawBline() { return DrawBline; }

  public:  // without description

    // inline void SetDrawLineWidth(G4double aVal)
    //   { TrajectoryVisAttributes.SetLineWidth(aVal); }
    // inline void SetDrawLineStyle(G4VisAttributes::LineStyle aStyle)
    //   { TrajectoryVisAttributes.SetLineStyle(aStyle); }
      // Future implementation...

  private:
  
   G4BlineTracer* fBlineTool;   
   G4Colour DrawColour;
   G4bool DrawBline;
   G4bool DrawPoints;
   G4double PointSize;
   std::vector<G4VisAttributes*> TrajectoryVisAttributes;
   std::vector<G4Polyline> TrajectoryPolyline;
   std::vector<G4Polymarker> TrajectoryPoints;
};

#endif
