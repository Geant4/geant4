// 
// G4Ellipsoid.hh,v 1.0 2005/02/25 12:00 gguerrie Exp
//
// class G4Ellipsoid
//
// A G4Ellipsoid is an ellipsoidal solid, optionally cut at two given z.
//
// Member Data:
//
//	xSemiAxis	semi-axis, x
//	ySemiAxis	semi-axis, y
//	zSemiAxis	semi-axis, z
//      zBottomCut      lower cut plane level, z (solid lies above this plane)
//      zTopCut         upper cut plane level, z (solid lies below this plane)
//
// First writing: G.Guerrieri, 25.02.2005
// Modified: G. Guerrieri, 08.06.2005

#ifndef G4Ellipsoid_HH
#define G4Ellipsoid_HH

#include "G4CSGSolid.hh"

class G4Ellipsoid : public G4CSGSolid {
public:
          G4Ellipsoid(const G4String& pName,
	           G4double pxSemiAxis,
	           G4double pySemiAxis,
	           G4double pzSemiAxis,
	           G4double pzBottomCut,
		   G4double pzTopCut);
		   
	  virtual ~G4Ellipsoid() ;
	  
	  // Access functions
		   
	  G4double    GetSemiAxisMax (int i) const  { return (i==0) ? xSemiAxis
						   : (i==1) ? ySemiAxis
						   : zSemiAxis; }
	  G4double    GetZBottomCut() const { return zBottomCut; }
	  G4double    GetZTopCut() const { return zTopCut; }

	  void   SetSemiAxis  (G4double newxSemiAxis,
			      G4double newySemiAxis,
			      G4double newzSemiAxis)
                   { xSemiAxis= newxSemiAxis; ySemiAxis= newySemiAxis; zSemiAxis= newzSemiAxis;
		     semiAxisMax = xSemiAxis > ySemiAxis ? xSemiAxis : ySemiAxis;
		     if (zSemiAxis > semiAxisMax) semiAxisMax= zSemiAxis;
		     if (zBottomCut < -zSemiAxis) zBottomCut = -zSemiAxis;
		     if (zTopCut > +zSemiAxis) zTopCut = +zSemiAxis;
		   }
          void   SetZCuts (G4double newzBottomCut, G4double newzTopCut)
                   { zBottomCut = newzBottomCut; zTopCut = newzTopCut;
		     if (zBottomCut < -zSemiAxis) zBottomCut = -zSemiAxis;
		     if (zTopCut > +zSemiAxis) zTopCut = +zSemiAxis;
		   }

          void ComputeDimensions(G4VPVParameterisation* p,
                                 const G4int n,
                                 const G4VPhysicalVolume* pRep);

          G4bool CalculateExtent(const EAxis pAxis,
			         const G4VoxelLimits& pVoxelLimit,
			         const G4AffineTransform& pTransform,
			         G4double& pmin, G4double& pmax) const;
				 
          EInside Inside(const G4ThreeVector& p) const;

          G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

          G4double DistanceToIn(const G4ThreeVector& p,
	                        const G4ThreeVector& v) const;
	  
          G4double DistanceToIn(const G4ThreeVector& p) const;
	  
          G4double DistanceToOut(const G4ThreeVector& p,
	                         const G4ThreeVector& v,
			         const G4bool calcNorm=G4bool(false),
			         G4bool *validNorm=0,
				 G4ThreeVector *n=0) const;
				 
          G4double DistanceToOut(const G4ThreeVector& p) const;

             // Naming method (pseudo-RTTI : run-time type identification)

          virtual G4GeometryType  GetEntityType() const
                    { return G4String("G4Ellipsoid"); }

              // Visualisation functions
  
          void   DescribeYourselfTo(G4VGraphicsScene& scene) const;
       
          G4VisExtent GetExtent() const;
       
          G4Polyhedron* CreatePolyhedron() const;
       
          G4NURBS*      CreateNURBS() const;
       
protected:
 
          G4ThreeVectorList*
          CreateRotatedVertices(const G4AffineTransform& pTransform,
	                        G4int& noPolygonVertices) const;
	
private:

    G4double xSemiAxis,
	     ySemiAxis,
	     zSemiAxis,
	     semiAxisMax,
	     zBottomCut,
 	     zTopCut;
};
   	
#endif
