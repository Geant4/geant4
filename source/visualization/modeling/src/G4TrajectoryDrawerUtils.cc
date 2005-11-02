#include "G4TrajectoryDrawerUtils.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"

namespace G4TrajectoryDrawerUtils {

  void GetPoints(const G4VTrajectory& traj, G4Polyline& trajectoryLine,
		 G4Polymarker& auxiliaryPoints, G4Polymarker& stepPoints) 
  {
    for (G4int i=0; i<traj.GetPointEntries(); i++) {
      G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(i);
      
      const std::vector<G4ThreeVector>* auxiliaries
	= aTrajectoryPoint->GetAuxiliaryPoints();
      
      if (auxiliaries) {
	for (size_t iAux=0; iAux<auxiliaries->size(); ++iAux) {
	  const G4ThreeVector pos((*auxiliaries)[iAux]);
	  trajectoryLine.push_back(pos);
	  auxiliaryPoints.push_back(pos);
	}
      }
      const G4ThreeVector pos(aTrajectoryPoint->GetPosition());
      trajectoryLine.push_back(pos);
      stepPoints.push_back(pos);
    }    
  }
  
}
