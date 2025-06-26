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
// Jane Tinslay, John Allison, Joseph Perl November 2005
//
#include "G4TrajectoryDrawerUtils.hh"

#include "G4AttValue.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4UIcommand.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4RichTrajectoryPoint.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4VisTrajContext.hh"
#include "G4VModel.hh"
#include "G4ModelingParameters.hh"

#include <utility>

#define G4warn G4cout

namespace G4TrajectoryDrawerUtils
{

enum TimesValidity
{
  InvalidTimes,
  ValidTimes
};

TimesValidity GetPointsAndTimes(const G4VTrajectory& traj, const G4VisTrajContext& context,
                                G4Polyline& trajectoryLine, G4Polymarker& auxiliaryPoints,
                                G4Polymarker& stepPoints,
                                std::vector<G4double>& trajectoryLineTimes,
                                std::vector<G4double>& auxiliaryPointTimes,
                                std::vector<G4double>& stepPointTimes)
{
  TimesValidity validity = InvalidTimes;
  if (context.GetTimeSliceInterval()) validity = ValidTimes;

  // Memory for last trajectory point position for auxiliary point
  // time interpolation algorithm.  There are no auxiliary points
  // for the first trajectory point, so its initial value is
  // immaterial.
  G4ThreeVector lastTrajectoryPointPosition;

  // Keep positions.  Don't store unless first or different.
  std::vector<G4ThreeVector> positions;

  for (G4int iPoint = 0; iPoint < traj.GetPointEntries(); iPoint++) {
    G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(iPoint);
    const G4ThreeVector& trajectoryPointPosition = aTrajectoryPoint->GetPosition();

    // Only store if first or if different
    if (positions.size() == 0 || trajectoryPointPosition != positions[positions.size() - 1]) {
      // Pre- and Post-Point times from the trajectory point...
      G4double trajectoryPointPreTime = -std::numeric_limits<double>::max();
      G4double trajectoryPointPostTime = std::numeric_limits<double>::max();

      if (context.GetTimeSliceInterval() && validity == ValidTimes) {

        const auto richTrajectoryPoint = dynamic_cast<G4RichTrajectoryPoint*>(aTrajectoryPoint);
        if (richTrajectoryPoint) {
          trajectoryPointPreTime  = richTrajectoryPoint->GetPreStepPointGlobalTime();
          trajectoryPointPostTime = richTrajectoryPoint->GetPostStepPointGlobalTime();
          if (trajectoryPointPostTime > 1000000. * context.GetTimeSliceInterval()) {
            static G4bool warnedTimeExceeded = false;
            if (!warnedTimeExceeded) {
              G4Exception("G4TrajectoryUtils::GetPointsAndTimes", "modeling0128", JustWarning,
                          "Trajectory point times exceed limit (>1000000 time slice)."
                          "\nNormal trajectory will be drawn.");
              warnedTimeExceeded = true;
            }
            validity = InvalidTimes;
          }
        } else {
          static G4bool warnedNoAttValues = false;
          if (!warnedNoAttValues) {
            G4Exception("G4TrajectoryUtils::GetPointsAndTimes", "modeling0127", JustWarning,
                        "Trajectory point times not found. Use G4RichTrajectory:"
                        "\n  /vis/scene/add/trajectories rich");
            warnedNoAttValues = true;
          }
          validity = InvalidTimes;
        }
      }

      const std::vector<G4ThreeVector>* auxiliaries = aTrajectoryPoint->GetAuxiliaryPoints();
      if (0 != auxiliaries) {
        for (size_t iAux = 0; iAux < auxiliaries->size(); ++iAux) {
          const G4ThreeVector& auxPointPosition = (*auxiliaries)[iAux];
          if (positions.size() == 0 || auxPointPosition != positions[positions.size() - 1]) {
            // Only store if first or if different
            positions.push_back(trajectoryPointPosition);
            trajectoryLine.push_back(auxPointPosition);
            auxiliaryPoints.push_back(auxPointPosition);
            if (validity == ValidTimes) {
              // Interpolate time for auxiliary points...
              G4double s1 = (auxPointPosition - lastTrajectoryPointPosition).mag();
              G4double s2 = (trajectoryPointPosition - auxPointPosition).mag();
              G4double t = trajectoryPointPreTime
                           + (trajectoryPointPostTime - trajectoryPointPreTime) * (s1 / (s1 + s2));
              trajectoryLineTimes.push_back(t);
              auxiliaryPointTimes.push_back(t);
            }
          }
        }
      }

      positions.push_back(trajectoryPointPosition);
      trajectoryLine.push_back(trajectoryPointPosition);
      stepPoints.push_back(trajectoryPointPosition);
      if (validity == ValidTimes) {
        trajectoryLineTimes.push_back(trajectoryPointPostTime);
        stepPointTimes.push_back(trajectoryPointPostTime);
      }
      lastTrajectoryPointPosition = trajectoryPointPosition;
    }
  }
  return validity;
}

static void SliceLine(G4double timeIncrement, G4Polyline& trajectoryLine,
                      std::vector<G4double>& trajectoryLineTimes)
{
  // Assumes valid arguments from GetPoints and GetTimes.

  G4Polyline newTrajectoryLine;
  std::vector<G4double> newTrajectoryLineTimes;

  newTrajectoryLine.push_back(trajectoryLine[0]);
  newTrajectoryLineTimes.push_back(trajectoryLineTimes[0]);
  size_t lineSize = trajectoryLine.size();
  if (lineSize > 1) {
    for (size_t i = 1; i < trajectoryLine.size(); ++i) {
      G4double deltaT = trajectoryLineTimes[i] - trajectoryLineTimes[i - 1];
      if (deltaT > 0.) {
        G4double practicalTimeIncrement = std::max(timeIncrement, deltaT / 100.);
        for (G4double t = (int(trajectoryLineTimes[i - 1] / practicalTimeIncrement) + 1)
                          * practicalTimeIncrement;
             t <= trajectoryLineTimes[i]; t += practicalTimeIncrement)
        {
          G4ThreeVector pos = trajectoryLine[i - 1]
                              + (trajectoryLine[i] - trajectoryLine[i - 1])
                                  * ((t - trajectoryLineTimes[i - 1]) / deltaT);
          newTrajectoryLine.push_back(pos);
          newTrajectoryLineTimes.push_back(t);
        }
      }
      newTrajectoryLine.push_back(trajectoryLine[i]);
      newTrajectoryLineTimes.push_back(trajectoryLineTimes[i]);
    }
  }

  trajectoryLine = std::move(newTrajectoryLine);
  trajectoryLineTimes = std::move(newTrajectoryLineTimes);
}

static void DrawWithoutTime(const G4VisTrajContext& myContext, G4Polyline& trajectoryLine,
                            G4Polymarker& auxiliaryPoints, G4Polymarker& stepPoints)
{
  // Draw without time slice information

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (0 == pVVisManager) return;

  if (myContext.GetDrawLine() && myContext.GetLineVisible()) {
    G4VisAttributes trajectoryLineAttribs(myContext.GetLineColour());
    trajectoryLineAttribs.SetLineWidth(myContext.GetLineWidth());
    trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);

    pVVisManager->Draw(trajectoryLine);
  }

  if (myContext.GetDrawAuxPts() && myContext.GetAuxPtsVisible() && (auxiliaryPoints.size() > 0)) {
    auxiliaryPoints.SetMarkerType(myContext.GetAuxPtsType());
    auxiliaryPoints.SetSize(myContext.GetAuxPtsSizeType(), myContext.GetAuxPtsSize());
    auxiliaryPoints.SetFillStyle(myContext.GetAuxPtsFillStyle());

    G4VisAttributes auxiliaryPointsAttribs(myContext.GetAuxPtsColour());
    auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);

    pVVisManager->Draw(auxiliaryPoints);
  }

  if (myContext.GetDrawStepPts() && myContext.GetStepPtsVisible() && (stepPoints.size() > 0)) {
    stepPoints.SetMarkerType(myContext.GetStepPtsType());
    stepPoints.SetSize(myContext.GetStepPtsSizeType(), myContext.GetStepPtsSize());
    stepPoints.SetFillStyle(myContext.GetStepPtsFillStyle());

    G4VisAttributes stepPointsAttribs(myContext.GetStepPtsColour());
    stepPoints.SetVisAttributes(&stepPointsAttribs);

    pVVisManager->Draw(stepPoints);
  }
}

static void DrawWithTime(const G4VisTrajContext& myContext, G4Polyline& trajectoryLine,
                         G4Polymarker& auxiliaryPoints, G4Polymarker& stepPoints,
                         std::vector<G4double>& trajectoryLineTimes,
                         std::vector<G4double>& auxiliaryPointTimes,
                         std::vector<G4double>& stepPointTimes)
{
  // Draw with time slice information

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (0 == pVVisManager) return;

  const auto& timeParameters = G4VModel::GetCurrentModelingParameters()->GetTimeParameters();
  const auto& viewerStartTime    = timeParameters.fStartTime;
  const auto& viewerEndTime      = timeParameters.fEndTime;
  const auto& viewerFadeFactor   = timeParameters.fFadeFactor;

  if (myContext.GetDrawLine() && myContext.GetLineVisible()) {
    G4VisAttributes trajectoryLineAttribs(myContext.GetLineColour());
    trajectoryLineAttribs.SetLineWidth(myContext.GetLineWidth());

    for (size_t i = 1; i < trajectoryLine.size(); ++i) {
      const auto& sliceStartTime = trajectoryLineTimes[i - 1];
      const auto& sliceEndTime = trajectoryLineTimes[i];
      if (sliceStartTime >= viewerStartTime && sliceEndTime <= viewerEndTime) {
        G4Polyline slice;
        slice.push_back(trajectoryLine[i - 1]);
        slice.push_back(trajectoryLine[i]);
        trajectoryLineAttribs.SetStartTime(trajectoryLineTimes[i - 1]);
        trajectoryLineAttribs.SetEndTime(trajectoryLineTimes[i]);
        if (viewerStartTime >= 0.) {
          auto opacityMultiplier = 1. - viewerFadeFactor
          * (viewerEndTime - trajectoryLineTimes[i]) / (viewerEndTime - viewerStartTime);
          auto colour = myContext.GetLineColour();
          colour.SetAlpha(opacityMultiplier * colour.GetAlpha());
          trajectoryLineAttribs.SetColour(colour);
        }
        slice.SetVisAttributes(&trajectoryLineAttribs);
        pVVisManager->Draw(slice);
      }
    }
  }

  if (myContext.GetDrawAuxPts() && myContext.GetAuxPtsVisible() && (auxiliaryPoints.size() > 0)) {
    G4VisAttributes auxiliaryPointsAttribs(myContext.GetAuxPtsColour());

    for (size_t i = 0; i < auxiliaryPoints.size(); ++i) {
      const auto& auxStartTime = auxiliaryPointTimes[i];
      const auto& auxEndTime = auxiliaryPointTimes[i];
      if (auxStartTime >= viewerStartTime && auxEndTime <= viewerEndTime) {
        G4Polymarker point;
        point.push_back(auxiliaryPoints[i]);
        point.SetMarkerType(myContext.GetAuxPtsType());
        point.SetSize(myContext.GetAuxPtsSizeType(), myContext.GetAuxPtsSize());
        point.SetFillStyle(myContext.GetAuxPtsFillStyle());
        auxiliaryPointsAttribs.SetStartTime(auxiliaryPointTimes[i]);
        auxiliaryPointsAttribs.SetEndTime(auxiliaryPointTimes[i]);
        if (viewerStartTime >= 0.) {
          auto opacityMultiplier = 1. - viewerFadeFactor
          * (viewerEndTime - auxiliaryPointTimes[i] ) / (viewerEndTime - viewerStartTime);
          auto colour = myContext.GetAuxPtsColour();
          colour.SetAlpha(opacityMultiplier * colour.GetAlpha());
          auxiliaryPointsAttribs.SetColour(colour);
        }
        point.SetVisAttributes(&auxiliaryPointsAttribs);
        pVVisManager->Draw(point);
      }
    }
  }

  if (myContext.GetDrawStepPts() && myContext.GetStepPtsVisible() && (stepPoints.size() > 0)) {
    G4VisAttributes stepPointsAttribs(myContext.GetStepPtsColour());

    for (size_t i = 0; i < stepPoints.size(); ++i) {
      const auto& pointStartTime = stepPointTimes[i];
      const auto& pointEndTime = stepPointTimes[i];
      if (pointStartTime >= viewerStartTime && pointEndTime <= viewerEndTime) {
        G4Polymarker point;
        point.push_back(stepPoints[i]);
        point.SetMarkerType(myContext.GetStepPtsType());
        point.SetSize(myContext.GetStepPtsSizeType(), myContext.GetStepPtsSize());
        point.SetFillStyle(myContext.GetStepPtsFillStyle());
        stepPointsAttribs.SetStartTime(stepPointTimes[i]);
        stepPointsAttribs.SetEndTime(stepPointTimes[i]);
        if (viewerStartTime >= 0.) {
          auto opacityMultiplier = 1. - viewerFadeFactor
          * (viewerEndTime - stepPointTimes[i]) / (viewerEndTime - viewerStartTime);
          auto colour = myContext.GetStepPtsColour();
          colour.SetAlpha(opacityMultiplier * colour.GetAlpha());
          stepPointsAttribs.SetColour(colour);
        }
        point.SetVisAttributes(&stepPointsAttribs);
        pVVisManager->Draw(point);
      }
    }
  }
}

void DrawLineAndPoints(const G4VTrajectory& traj, const G4VisTrajContext& context)
{
  // Return if don't need to do anything
  if (!context.GetDrawLine() && !context.GetDrawAuxPts() && !context.GetDrawStepPts()) return;

  // Get points and times (times are returned only if time-slicing
  // is requested).
  G4Polyline trajectoryLine;
  G4Polymarker stepPoints;
  G4Polymarker auxiliaryPoints;
  std::vector<G4double> trajectoryLineTimes;
  std::vector<G4double> stepPointTimes;
  std::vector<G4double> auxiliaryPointTimes;

  TimesValidity validity =
    GetPointsAndTimes(traj, context, trajectoryLine, auxiliaryPoints, stepPoints,
                      trajectoryLineTimes, auxiliaryPointTimes, stepPointTimes);

  if (validity == ValidTimes) {
    SliceLine(context.GetTimeSliceInterval(), trajectoryLine, trajectoryLineTimes);

    DrawWithTime(context, trajectoryLine, auxiliaryPoints, stepPoints, trajectoryLineTimes,
                 auxiliaryPointTimes, stepPointTimes);
  }
  else {
    DrawWithoutTime(context, trajectoryLine, auxiliaryPoints, stepPoints);
  }
}
}  // namespace G4TrajectoryDrawerUtils
