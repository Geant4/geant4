//
// Demonstrate some issues with G4Tubs
//

#include <assert.h>
#include <iomanip>
#include <math.h>

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"

#include "G4Tubs.hh"

// #include "G4Tubs13.hh"

#include "G4Hype.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"


const G4String OutputInside(const EInside a)
{
  switch(a) 
  {
    case kInside:  return " Inside"; 
    case kOutside: return "Outside";
    case kSurface: return "Surface";
  }
  return "   ????";
}


int main()
{

	
// Voxel problems
//
// Make a simple cylinder
	
  G4Tubs test3( "simple", 0.5*m, 1*m, 0.1*m, 0, 2*pi );
	
	
// Make an equivalent hyperbolic tube
	
  G4Hype test4( "simple2", 0.5*m, 1*m, 0.0, 0.0, 0.1*m );
	
	
// Build a transform that is a simple rotation
// of 85 degrees around the z axis
	
  G4RotationMatrix rotate;

  //  rotate.rotateZ(85.0*degree);
  //  rotate.rotateZ(90.0*degree);

  rotate.rotateZ(270.0*degree);
  G4AffineTransform transform( rotate );
	
	
// Here is a voxel limit that extends into the
// middle of the solid
	
  G4VoxelLimits voxel1;
  voxel1.AddLimit( kXAxis, -0.1*m, 0.1*m );
  voxel1.AddLimit( kYAxis,  0.8*m, 2.0*m );
  voxel1.AddLimit( kZAxis, -0.001*m, 0.001*m );
	
	
// Test 1 CalculateExtent. Answer should be:
//          Min = 0.8
//          Max = something just larger than 1*m
	
  G4bool answer3, answer4;
  G4double min3, max3, min4, max4;
	
  answer3 = test3.CalculateExtent( kYAxis, voxel1, transform, min3, max3 );
  answer4 = test4.CalculateExtent( kYAxis, voxel1, transform, min4, max4 );
	
  G4cout << "-----------------" << G4std::endl;
  G4cout << "Voxel test 1 " << G4std::endl;
  G4cout << "Should be: "
	 << G4std::setw(7) << true
	 << G4std::setw(10) << 800
	 << G4std::setw(10) << ">1000" << G4std::endl;
  G4cout << "G4Tubs:    "
	 << G4std::setw(7) << answer3
	 << G4std::setw(10) << min3
	 << G4std::setw(10) << max3 << G4std::endl;
  G4cout << "G4Hype:    "
	 << G4std::setw(7) << answer4
	 << G4std::setw(10) << min4
	 << G4std::setw(10) << max4 << G4std::endl;
	
	
// Here is a voxel limit that starts and ends in the solid
	
  G4VoxelLimits voxel2;
  // voxel2.AddLimit( kXAxis, -0.1*m, 0.1*m );
  voxel2.AddLimit( kXAxis, -0.8*m, 0.8*m );
  voxel2.AddLimit( kYAxis, -0.8*m, 0.8*m );
  voxel2.AddLimit( kZAxis, -0.001*m, 0.001*m );
	
	
// Test 2 CalculateExtent. Answer should be:
//          Min = -0.8
//          Max = +0.8
	
//  answer3 = test3.CalculateExtent( kXAxis, voxel2, transform, min3, max3 );
//  answer4 = test4.CalculateExtent( kXAxis, voxel2, transform, min4, max4 );

  answer3 = test3.CalculateExtent( kYAxis, voxel2, transform, min3, max3 );
  answer4 = test4.CalculateExtent( kYAxis, voxel2, transform, min4, max4 );
	
	G4cout << "-----------------" << G4std::endl;
	G4cout << "Voxel test 2 " << G4std::endl;
	G4cout << "Should be: "
	       << G4std::setw(7) << true
	       << G4std::setw(10) << -800
	       << G4std::setw(10) << 800 << G4std::endl;
	G4cout << "G4Tubs:    "
	       << G4std::setw(7) << answer3
	       << G4std::setw(10) << min3
	       << G4std::setw(10) << max3 << G4std::endl;
	G4cout << "G4Hype:    "
	       << G4std::setw(7) << answer4
	       << G4std::setw(10) << min4
	       << G4std::setw(10) << max4 << G4std::endl;
	
	
// Here is a tiny voxel that just skims the solid.
	
  G4VoxelLimits voxel3;
  // voxel3.AddLimit( kXAxis, -0.001*m,  0.001*m );

  voxel3.AddLimit( kXAxis, 0.998*m,  1.003*m );
  voxel3.AddLimit( kYAxis,  0.998*m,  1.003*m );
  voxel3.AddLimit( kZAxis, -0.001*m,  0.001*m );
	
	
// Test 3 CalculateExtent. Answer should be:
//           Min = 0.998
//           Max = Something larger than 1.0
//

// answer3 = test3.CalculateExtent( kXAxis, voxel3, transform, min3, max3 );
// answer4 = test4.CalculateExtent( kXAxis, voxel3, transform, min4, max4 );
  answer3 = test3.CalculateExtent( kYAxis, voxel3, transform, min3, max3 );
  answer4 = test4.CalculateExtent( kYAxis, voxel3, transform, min4, max4 );
	
        G4cout << "-----------------" << G4std::endl;
	G4cout << "Voxel test 3 " << G4std::endl;
	G4cout << "Should be: "
	       << G4std::setw(7) << true
	       << G4std::setw(10) << 998
	       << G4std::setw(10) << ">1000" << G4std::endl;
	G4cout << "G4Tubs:    "
	       << G4std::setw(7) << answer3
	       << G4std::setw(10) << min3
	       << G4std::setw(10) << max3 << G4std::endl;
	G4cout << "G4Hype:    "
	       << G4std::setw(7) << answer4
	       << G4std::setw(10) << min4
	       << G4std::setw(10) << max4 << G4std::endl;
}
