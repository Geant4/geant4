#ifndef G4TRAJECTORYDRAWERUTILS
#define G4TRAJECTORYDRAWERUTILS

class G4Polyline;
class G4Polymarker;
class G4VTrajectory;

namespace G4TrajectoryDrawerUtils {

  void GetPoints(const G4VTrajectory& traj, G4Polyline& trajectoryLine,
		 G4Polymarker& auxiliaryPoints, G4Polymarker& stepPoints);

}

#endif
