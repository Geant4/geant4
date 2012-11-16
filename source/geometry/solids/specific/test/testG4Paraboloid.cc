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
#undef NDEBUG
#include <assert.h>

#include "G4GeometryTolerance.hh"
#include "G4Paraboloid.hh"
#include "G4Polyhedron.hh"
#include "Randomize.hh"
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
  // G4cout<<paraboloid1.GetPolyhedron()->GetSurfaceArea()<<G4endl;
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
  

  // Add Test Distance To Out for User Problem with Substraction Of Paraboloid :: Problem Report N1015
   G4Paraboloid paraboloid3("Paraboloid1", 72., 0., 240.);
  //  G4Paraboloid paraboloid3("Paraboloid1", 72., 210., 240.);//This test is Ok, Intersection with DZ working correctly
   G4double distOut1;
   distOut1=paraboloid3.DistanceToOut(G4ThreeVector(200., 0., 72.), G4ThreeVector(0.0,0.,-1.), false, NULL, NULL);
   assert((distOut1-44.0)<1e-5);
   //G4cout<<"Test3::distOut="<<distOut1<<G4endl;

   distOut1=paraboloid3.DistanceToOut(G4ThreeVector(200., 0., 72.), G4ThreeVector(0.000000000001,0.,-1.), false, NULL, NULL);
    assert((distOut1-44.0)<1e-5);
   //G4cout<<"Test3::distOut="<<distOut1<<G4endl;
   // Add Rotation On the Surface around "Problem" Point(200.,0.,72.) 

    G4double tolA=0.25e-7; 
    G4ThreeVector In (200.,0.,72.);
    EInside in; 
    G4ThreeVector Dir,vDir,pSurf;
    Dir=G4ThreeVector(10*tolA,0.,-1.);
    vDir=Dir.unit(); 
    in=paraboloid3.Inside(In);
    if(in!=kSurface)G4cout<<"Problem ::Point p"<<In<<" is NOT On Surface"<<G4endl;
    
    for(G4int i=0; i<100; i++){
     
      G4double distOut=paraboloid3.DistanceToOut(In,vDir);
      G4cout.precision(16);
      pSurf=In+vDir*distOut;
      in=paraboloid3.Inside(pSurf);
      if(in!=kSurface) G4cout<<"Problem :: ps is Not on the Surface DistToOut="<<distOut<<" dir="<<vDir<<G4endl;
      Dir=Dir-G4ThreeVector(tolA,0,0);
      vDir=Dir.unit();

   }
    // Rotation from the Lower DZ plane
     Dir=G4ThreeVector(10*tolA,0.,1.);
     vDir=Dir.unit(); 
     In= G4ThreeVector (0.,0.,-72.);
     for(G4int i=0; i<100; i++){
     
      G4double distOut=paraboloid3.DistanceToOut(In,vDir);
      G4cout.precision(16);
      pSurf=In+vDir*distOut;
      in=paraboloid3.Inside(pSurf);
      if(in!=kSurface) G4cout<<"Problem :: ps is Not on the Surface DistToOut="<<distOut<<" dir="<<vDir<<G4endl;
      Dir=Dir-G4ThreeVector(tolA,0,0);
      vDir=Dir.unit();
     }
  // Paraboloid with R1!=0 and with Point situated On both Surfaces : Z and Parabolic
   G4Paraboloid paraboloid5("Paraboloid1", 10., 0., 100.);
   distOut1=paraboloid5.DistanceToOut(G4ThreeVector(100., 0., 10.), G4ThreeVector(0.0,0.,-1.), false, NULL, NULL);
   assert((distOut1-0.0)<1e-5);
  
    Dir=G4ThreeVector(-0.9,0.,-0.1);
    vDir=Dir.unit();
    distOut1=paraboloid5.DistanceToOut(G4ThreeVector(100., 0., 10.), vDir, false, NULL, NULL);
    pSurf=G4ThreeVector(100.0,0.,10.)+vDir*distOut1;
    in=paraboloid5.Inside(pSurf);
    if(in!=kSurface)G4cout<<"Problem :: ps is Not on the Surface "<<pSurf<<" vDir="<<vDir<<G4endl;

  return 0;
}
