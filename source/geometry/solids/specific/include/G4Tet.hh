//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of     *
// * the Vanderbilt University Free Electron Laser Center             *
// * Vanderbilt University, Nashville, TN, USA                        *
// * development supported by:                                        *
// * United States MFEL program  under grant FA9550-04-1-0045         *
// * and NASA under contract number NNG04CT05P                        *
// * and was written by Marcus H. Mendenhall and Robert A. Weller     *
// *                                                                  *
// * contributed to the Geant4 Core, January, 2005                    *
// *                                                                  *
// ********************************************************************

#ifndef G4TET_HH
#define G4TET_HH

//JA #include "G4CSGSolid.hh"
#include "G4VSolid.hh"

class G4Tet : public G4VSolid   //JA  was G4CSGSolid 
{

public:  // with description

       G4Tet(const G4String& pName, 
	           G4ThreeVector anchor,
	           G4ThreeVector p2,
		   G4ThreeVector p3,
		   G4ThreeVector p4, 
	     G4bool *degeneracyFlag=0);

        virtual ~G4Tet();

public:
	static G4bool CheckDegeneracy(
				G4ThreeVector anchor,
				G4ThreeVector p2,
				G4ThreeVector p3,
				G4ThreeVector p4);

    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const;

	void PrintWarnings(G4bool flag) { warningFlag=flag; };
	
  // Accessors and modifiers

  // Methods for solid

    EInside Inside(const G4ThreeVector& p) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const;

    G4double DistanceToIn(const G4ThreeVector& p) const;

    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm=false,
                                 G4bool *validNorm=0, G4ThreeVector *n=0) const;

    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4ThreeVector GetPointOnFace(G4ThreeVector p1, G4ThreeVector p2, 
				 G4ThreeVector p3, G4double& area) const;

    G4ThreeVector GetPointOnSurface() const;

    G4GeometryType GetEntityType() const;

    std::ostream& StreamInfo(std::ostream& os) const;

  // Functions for visualization

    void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent   GetExtent          () const;
    G4Polyhedron* CreatePolyhedron   () const;
    G4NURBS*      CreateNURBS        () const;

  private:
	static const char CVSVers[];

  public:   // without description 
	const char *CVSHeaderVers() { return 
	   "$Id: G4Tet.hh,v 1.2 2005-08-03 15:53:42 danninos Exp $";
	}

	const char *CVSFileVers() { return CVSVers; }

  protected:  // with description

    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;
      // Create the List of transformed vertices in the format required
      // for G4VSolid:: ClipCrossSection and ClipBetweenSections.

  private:

    G4ThreeVector fAnchor, fP2, fP3, fP4, fMiddle;
	G4ThreeVector fNormal123, fNormal142, fNormal134, fNormal234;
	G4bool warningFlag;
	
	G4double fCdotN123, fCdotN142, fCdotN134, fCdotN234;
	
	G4double fXMin, fXMax, fYMin, fYMax, fZMin, fZMax;
	G4double fDx, fDy, fDz, fTol, fMaxSize;
};

#endif
