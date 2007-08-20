#include <assert.h>

#include "G4GeometryTolerance.hh"
#include "G4Paraboloid.hh"
#include "G4Polyhedron.hh"
int main()
{
  G4ThreeVector testPoint;

  G4Paraboloid paraboloid1("Paraboloid1", 1., 0., 3.);
  G4Paraboloid paraboloid2("Paraboloid2", .01, 2., 3.);

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // Test Inside function.
  assert(paraboloid1.Inside(G4ThreeVector(2.598076211, 1.5, 1.)) == kSurface);
  assert(paraboloid1.Inside(G4ThreeVector((2 * kCarTolerance + 3.) * std::sqrt(3.) / 2., (2. * kCarTolerance + 3.) / 2., 1.)) == kOutside);
  assert(paraboloid1.Inside(G4ThreeVector((- 2 * kCarTolerance + 3.) * std::sqrt(3.) / 2., (- 2. * kCarTolerance + 3.) / 2., 1.)) == kSurface);

  assert(paraboloid1.Inside(G4ThreeVector(std::sqrt(4.5) / std::sqrt(2.), std::sqrt(4.5) / std::sqrt(2.), 0.)) == kSurface);
  assert(paraboloid1.Inside(G4ThreeVector((  2 * kCarTolerance + std::sqrt(4.5)) / std::sqrt(2.), (  2. * kCarTolerance + std::sqrt(4.5)) / std::sqrt(2.), 0.)) == kOutside);
  assert(paraboloid1.Inside(G4ThreeVector((- 2 * kCarTolerance + std::sqrt(4.5)) / std::sqrt(2.), (- 2. * kCarTolerance + std::sqrt(4.5)) / std::sqrt(2.), 0.)) == kInside);

  assert(paraboloid1.Inside(G4ThreeVector(0., 0., -1.)) == kSurface);
  assert(paraboloid1.Inside(G4ThreeVector(0., 0., 1.)) == kSurface);
  assert(paraboloid1.Inside(G4ThreeVector(0., 0., - 1. - 2 * kCarTolerance)) == kOutside);
  assert(paraboloid1.Inside(G4ThreeVector(0., 0.,   1. + 2 * kCarTolerance)) == kOutside);
  assert(paraboloid1.Inside(G4ThreeVector(0., 0., - 1. + 2 * kCarTolerance)) == kInside);
  assert(paraboloid1.Inside(G4ThreeVector(0., 0.,   1. - 2 * kCarTolerance)) == kInside);
  assert(paraboloid1.Inside(G4ThreeVector((- 2 * kCarTolerance) / std::sqrt(2.), ( 2. * kCarTolerance ) / std::sqrt(2.), -1.)) == kOutside);

  assert(paraboloid2.Inside(G4ThreeVector(2.598076211, 1.5, .01)) == kSurface);
  assert(paraboloid2.Inside(G4ThreeVector((2 * kCarTolerance + 3.) * std::sqrt(3.) / 2., (2. * kCarTolerance + 3.) / 2., .01)) == kOutside);
  assert(paraboloid2.Inside(G4ThreeVector((- 2 * kCarTolerance + 3.) * std::sqrt(3.) / 2., (- 2. * kCarTolerance + 3.) / 2., .01)) == kSurface);

  assert(paraboloid2.Inside(G4ThreeVector(std::sqrt(6.5) / std::sqrt(2.), std::sqrt(6.5) / std::sqrt(2.), 0.)) == kSurface);
  assert(paraboloid2.Inside(G4ThreeVector((  2. * kCarTolerance + std::sqrt(6.5)) / std::sqrt(2.), (  2. * kCarTolerance + std::sqrt(6.5)) / std::sqrt(2.), 0.)) == kOutside);
  assert(paraboloid2.Inside(G4ThreeVector((- 2. * kCarTolerance + std::sqrt(6.5)) / std::sqrt(2.), (- 2. * kCarTolerance + std::sqrt(6.5)) / std::sqrt(2.), 0.)) == kInside);

  assert(paraboloid2.Inside(G4ThreeVector(0., 0., -.01)) == kSurface);
  assert(paraboloid2.Inside(G4ThreeVector(0., 0., .01)) == kSurface);
  assert(paraboloid2.Inside(G4ThreeVector(0., 0., - .01 - 2 * kCarTolerance)) == kOutside);
  assert(paraboloid2.Inside(G4ThreeVector(0., 0.,   .01 + 2 * kCarTolerance)) == kOutside);
  assert(paraboloid2.Inside(G4ThreeVector(0., 0., - .01 + 2 * kCarTolerance)) == kInside);
  assert(paraboloid2.Inside(G4ThreeVector(0., 0.,   .01 - 2 * kCarTolerance)) == kInside);
  assert(paraboloid2.Inside(G4ThreeVector((- 2. - 2 * kCarTolerance) / std::sqrt(2.), (2. + 2. * kCarTolerance ) / std::sqrt(2.), -.01)) == kOutside);
  assert(paraboloid2.Inside(G4ThreeVector((- 2. + 2 * kCarTolerance) / std::sqrt(2.), (2. - 2. * kCarTolerance ) / std::sqrt(2.), -.01)) == kSurface);
  assert(paraboloid2.Inside(G4ThreeVector((- 2.) / std::sqrt(2.), (2.) / std::sqrt(2.), -.01)) == kSurface);

  // Test DistanceToIn(p, v) function.
  assert(std::fabs(paraboloid1.DistanceToIn( G4ThreeVector(3.,2.,10.), 
           G4ThreeVector(-0.099014754299999993558678568206233, -0.099014754299999993558678568206233, -0.99014754279999994679428709787317)) 
          - 9.0895544461477406628091557649896) < kCarTolerance);
  assert(std::fabs(paraboloid1.DistanceToIn(G4ThreeVector(1., 0., 1.), G4ThreeVector(0., 0., -1.))) < kCarTolerance);
  assert(std::fabs(paraboloid1.DistanceToIn(G4ThreeVector(0., 0., -2.), G4ThreeVector(0., 0., 1.)) - 1) < kCarTolerance);
  assert(std::fabs(paraboloid1.DistanceToIn(G4ThreeVector(-2., 0., -1.), G4ThreeVector(1/std::sqrt(2.), 0., 1/std::sqrt(2.))) - 1/ std::sqrt(2.)) < kCarTolerance);
  assert(std::fabs(paraboloid1.DistanceToIn(G4ThreeVector(1., 0., 1.), G4ThreeVector(1., 0., 0.)) ) >= kInfinity);
  assert(std::fabs(paraboloid1.DistanceToIn(G4ThreeVector(1., 0., 1.), G4ThreeVector(1., 0., -0.0000001)) ) <= kCarTolerance);

  assert(std::fabs(paraboloid2.DistanceToIn( G4ThreeVector(3.,2.,10.), 
           G4ThreeVector(-0.099014754299999993558678568206233, -0.099014754299999993558678568206233, -0.99014754279999994679428709787317)) 
          - 10.0894054352) < kCarTolerance);
  assert(std::fabs(paraboloid2.DistanceToIn(G4ThreeVector(1., 0., .01), G4ThreeVector(0., 0., -1.))) < kCarTolerance);
  assert(std::fabs(paraboloid2.DistanceToIn(G4ThreeVector(0., 0., -1.), G4ThreeVector(0., 0., 1.)) - 0.99) < kCarTolerance);
  assert(std::fabs(paraboloid2.DistanceToIn(G4ThreeVector(-2.2, 0., -.01), G4ThreeVector(1/std::sqrt(2.), 0., 1/std::sqrt(2.))) - 0.004669633692) < kCarTolerance);
  assert(std::fabs(paraboloid2.DistanceToIn(G4ThreeVector(1., 0., .01), G4ThreeVector(1., 0., 0.)) ) >= kInfinity);
  assert(std::fabs(paraboloid2.DistanceToIn(G4ThreeVector(1., 0., .01), G4ThreeVector(1., 0., -0.0000001)) ) <= kCarTolerance);

  // Test DistanceToOut(p, v, ...) function.
  assert(std::fabs(paraboloid1.DistanceToOut(G4ThreeVector(1., 0., 1.), G4ThreeVector(0., 0., -1.), false, NULL, NULL) - 1.7777777777777776) < kCarTolerance);
  assert(std::fabs(paraboloid1.DistanceToOut(G4ThreeVector(0., 0., 0.), G4ThreeVector(0., 0., -1.), false, NULL, NULL) - 1) < kCarTolerance);
  assert(std::fabs(paraboloid1.DistanceToOut(G4ThreeVector(0., 0., 0.), G4ThreeVector(1. / std::sqrt(2.), 1. / std::sqrt(2.), 0.), false, NULL, NULL) - 2.1213203435) < kCarTolerance);
  assert(std::fabs(paraboloid1.DistanceToOut(G4ThreeVector(0., 0., 0.), G4ThreeVector(1. / std::sqrt(2.), 0., 1. / std::sqrt(2.)), false, NULL, NULL) - 1.4142135623) < kCarTolerance);
  assert(std::fabs(paraboloid1.DistanceToOut(G4ThreeVector(0., 0., 0.), G4ThreeVector(1. / std::sqrt(2.), 0., -1. / std::sqrt(2.)), false, NULL, NULL) - 1.1912334058) < kCarTolerance);
  assert(std::fabs(paraboloid1.DistanceToOut(G4ThreeVector(1., 0., 1.), G4ThreeVector(1., 0., 0.), false, NULL, NULL) - 2.) < kCarTolerance);
  assert(std::fabs(paraboloid1.DistanceToOut(G4ThreeVector(1., 0., 1.), G4ThreeVector(1., 0., 0.00000001), false, NULL, NULL)) < kCarTolerance);

  assert(std::fabs(paraboloid2.DistanceToOut(G4ThreeVector(2.5, 0., .01), G4ThreeVector(0., 0., -1.), false, NULL, NULL) - 0.010999999999999999) < kCarTolerance);
  assert(std::fabs(paraboloid2.DistanceToOut(G4ThreeVector(0., 0., 0.), G4ThreeVector(0., 0., -1.), false, NULL, NULL) - 0.01) < kCarTolerance);
  assert(std::fabs(paraboloid2.DistanceToOut(G4ThreeVector(0., 0., 0.), G4ThreeVector(1. / std::sqrt(2.), 1. / std::sqrt(2.), 0.), false, NULL, NULL) - 2.5495097567) < kCarTolerance);
  assert(std::fabs(paraboloid2.DistanceToOut(G4ThreeVector(0., 0., 0.), G4ThreeVector(1. / std::sqrt(2.), 0., 1. / std::sqrt(2.)), false, NULL, NULL) - 0.01414213562) < kCarTolerance);
  assert(std::fabs(paraboloid2.DistanceToOut(G4ThreeVector(0., 0., 0.), G4ThreeVector(1. / std::sqrt(2.), 0., -1. / std::sqrt(2.)), false, NULL, NULL) - 0.01414213562) < kCarTolerance);
  assert(std::fabs(paraboloid2.DistanceToOut(G4ThreeVector(1., 0., -.01), G4ThreeVector(1., 0., 0.), false, NULL, NULL) - 1.) < kCarTolerance);
  assert(std::fabs(paraboloid2.DistanceToOut(G4ThreeVector(1., 0., .01), G4ThreeVector(1., 0., 0.00000001), false, NULL, NULL)) < kCarTolerance);

  // Test volume function.
  assert(std::fabs(paraboloid1.GetCubicVolume() - 28.274334) < 1e-6);
  assert(std::fabs(paraboloid2.GetCubicVolume() - 0.408407) < 1e-6);

  // Test area function.
  //G4cout<<paraboloid1.GetPolyhedron()->GetSurfaceArea()<<G4endl;
  assert(std::fabs(paraboloid1.GetSurfaceArea() - 66.758844) < 1e-6);
  assert(std::fabs(paraboloid2.GetSurfaceArea() - 56.551935) < 1e-6);

  // Test SurfaceNormal(p) function.
  assert((paraboloid1.SurfaceNormal(G4ThreeVector(1., 2., 1.)) - G4ThreeVector(0., 0., 1.)).r() < kCarTolerance);
  assert((paraboloid1.SurfaceNormal(G4ThreeVector(1.732050808, 2.449489743, 1.)) - G4ThreeVector(0.40824829,0.57735027,0.70710678)).r() < 1e-5);
  assert((paraboloid1.SurfaceNormal(G4ThreeVector(1.573213272, 1.573213272, 0.1)) - G4ThreeVector(0.49718308,0.49718308,-0.71106819)).r() < 1e-5);
  assert((paraboloid1.SurfaceNormal(G4ThreeVector(0., 0., -1.)) - G4ThreeVector(0., 0., -1.)).r() < kCarTolerance);

  assert((paraboloid2.SurfaceNormal(G4ThreeVector(1., 2., 1.)) - G4ThreeVector(0., 0., 1.)).r() < kCarTolerance);
  assert((paraboloid2.SurfaceNormal(G4ThreeVector(1.732050808, 2.449489743, 0.01)) - G4ThreeVector(0.40824829,0.57735027,0.70710678)).r() < 1e-5);
  assert((paraboloid2.SurfaceNormal(G4ThreeVector(1.658312395, 2., 0.001)) - G4ThreeVector(0.013263635,0.015996545,-0.99978407)).r() < 1e-5);
  assert((paraboloid2.SurfaceNormal(G4ThreeVector(0., 0., -.01)) - G4ThreeVector(0., 0., -1.)).r() < kCarTolerance);

  return 0;
}
