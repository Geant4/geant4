// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4AffineTransform.cc,v 1.3 1999-12-15 14:50:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include <assert.h>
#include "G4AffineTransform.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "ApproxEqual.hh"


G4bool testG4AffineTransform()
{
	G4ThreeVector zeroVec,xVec(1,0,0),xyzVec(1,1,1),xyzrotVec(-1,1,1);
	G4RotationMatrix identity,xRot;
// NOTE: xRot = rotation such that x axis->y axis & y axis->-x axis
	xRot.rotateZ(-M_PI*0.5);

	G4AffineTransform origin;
	assert(origin.NetRotation()==identity);
	assert(origin.NetTranslation()==zeroVec);
	assert(!origin.IsRotated());
	assert(!origin.IsTranslated());

	G4AffineTransform rotTf(xRot);
	assert(rotTf.NetRotation()==xRot);
	assert(rotTf.NetTranslation()==zeroVec);
	assert(rotTf.IsRotated());
	assert(!rotTf.IsTranslated());

	G4AffineTransform txTf(xyzVec);
	assert(txTf.NetRotation()==identity);
	assert(txTf.NetTranslation()==xyzVec);
	assert(!txTf.IsRotated());
	assert(txTf.IsTranslated());

	G4AffineTransform rtTf(xRot,xyzVec);
	assert(rtTf.NetRotation()==xRot);
	assert(rtTf.NetTranslation()==xyzVec);
	assert(rtTf.IsRotated());
	assert(rtTf.IsTranslated());

	G4AffineTransform copyTf(rtTf);
	assert(copyTf==rtTf);

	G4AffineTransform compoundTf1=rotTf*txTf;
	assert(compoundTf1==rtTf);

	G4AffineTransform compoundTf2(rotTf);
	compoundTf2*=txTf;
	assert(compoundTf2==rtTf);

	G4AffineTransform compoundTf3;
	compoundTf3.Product(rotTf,txTf);
	assert(compoundTf3==rtTf);

	G4AffineTransform compoundTf4;
	compoundTf4.InverseProduct(rtTf,txTf);
	assert(ApproxEqual(compoundTf4,rotTf));
	compoundTf4.InverseProduct(rtTf,rtTf);
	assert(ApproxEqual(compoundTf4,identity));

	G4AffineTransform compoundTf5;
	compoundTf5.Product(rotTf,rtTf);
	G4AffineTransform compoundTf6;
	compoundTf6.InverseProduct(compoundTf5,rtTf);
	assert(ApproxEqual(compoundTf6,rotTf));

	assert(ApproxEqual(rotTf.TransformPoint(xyzVec),xyzrotVec));
	assert(ApproxEqual(rotTf.TransformAxis(xyzVec),xyzrotVec));
	assert(ApproxEqual(txTf.TransformPoint(xyzVec),G4ThreeVector(2,2,2)));
	assert(txTf.TransformAxis(xyzVec)==xyzVec);
	assert(ApproxEqual(rtTf.TransformPoint(xVec),G4ThreeVector(1,2,1)));
	assert(ApproxEqual(rtTf.TransformAxis(xVec),G4ThreeVector(0,1,0)));

	G4ThreeVector vec(0,0,1);
	rtTf.ApplyPointTransform(vec);
	assert(ApproxEqual(vec,G4ThreeVector(1,1,2)));
	vec=G4ThreeVector(-1,2,-3);
	rtTf.ApplyAxisTransform(vec);
	assert(ApproxEqual(vec,G4ThreeVector(-2,-1,-3)));
	rtTf.ApplyPointTransform(vec);
	assert(ApproxEqual(vec,G4ThreeVector(2,-1,-2)));

	G4AffineTransform invTf=rtTf.Inverse();

#if 0
	G4ThreeVector forwV= rtTf.TransformPoint(xyzVec);
	G4ThreeVector backV= invTf.TransformPoint(forwV);
        
        G4ThreeVector diffV= xyzVec - backV; 
        G4cout << " Diff of xyzVec and backV is " << diffV << G4endl;
#endif
	assert(ApproxEqual(invTf.TransformPoint(rtTf.TransformPoint(xyzVec)),
			   xyzVec));

	invTf*=rtTf;
// Might need tolerant checking:
	assert(ApproxEqual(invTf,origin));

	invTf=rtTf;
	invTf.Invert();

        G4double MaxAbsDiff(const G4AffineTransform &tf1,
	                    const G4AffineTransform &tf2); 
        G4double maxabsdiff= MaxAbsDiff( invTf, rtTf.Inverse() );
        G4cout << "Max difference is " << maxabsdiff << G4endl;
	assert( maxabsdiff <= 1.e-12 );

        G4AffineTransform  rtTf_inv=rtTf.Inverse();
	assert(MaxAbsDiff( invTf, rtTf.Inverse()) <= 1.e-12 );
#if 0
	assert(invTf==rtTf_inv);
         
	assert(invTf==rtTf.Inverse());
#endif

	G4AffineTransform txTf2(xyzVec);
	txTf2+=xyzVec;
	assert(txTf2.NetRotation()==identity);
	assert(ApproxEqual(txTf2.NetTranslation(),xyzVec*2));
	assert(txTf2!=txTf);
	txTf2-=xyzVec;
	assert(txTf2.NetRotation()==identity);
	assert(ApproxEqual(txTf2.NetTranslation(),xyzVec));
	assert(ApproxEqual(txTf2,txTf));
	txTf2.SetNetRotation(xRot);
	assert(txTf2.NetRotation()==xRot);
	txTf2.SetNetTranslation(xyzVec*3);
	assert(txTf2.NetTranslation()==xyzVec*3);

	return true;
}

int main()
{

#ifdef NDEBUG
	G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
	assert(testG4AffineTransform());
	return 0;
}

G4double MaxAbsDiff(const G4AffineTransform &tf1,
	            const G4AffineTransform &tf2)
{
        G4double maxabs= 0.0;
	for (G4int i=0;i<15;i++)
		{
                  G4double absdiff; 

		  absdiff= fabs(tf1[i]-tf2[i]);
                  maxabs= G4std::max(absdiff, maxabs); 
		}
	return maxabs;
}
