// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4makevol.cc,v 1.1 1999-01-07 16:06:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G3toG4.hh"
#include "G3VolTable.hh"
#include "G3MedTable.hh"
#include "G4BREPSolidPCone.hh"
#include "G4BREPSolidPolyhedra.hh"
#include "G4Para.hh"

G4double G3Bound(const G4String& s, const G4double& low,
                   const G4double& high, G4double val)
{
    if (val<low){
        G4cout << "Warning (" << s << ") " << val << " reset to " << low << endl;
        val=low;
    }
    if (val>high){
        G4cout << "Warning (" << s << ") " << val << " reset to " << high << endl;
        val=high;
    }
    return val;
}

        
G4LogicalVolume* G4makevol(G4String vname, G4String shape, G4int nmed,
                           G4double Rpar[], G4int npar)
{
//    G4cout << "shape " << shape << " name " << vname << " nmed " << nmed
//         << " npar " << npar << " parameters: ";
//    {
//        for (G4int i=0; i<npar; i++){
//            G4cout << Rpar[i] << " ";
//        }
//        G4cout << endl;
//    }
    
    G4double rangehi[3];
    G4double rangelo[3];
    EAxis axis = kXAxis;
  
        // Create the solid
    G4VSolid *solid = NULL;
  
    if ( shape == "BOX" ) {
        G4double pX = Rpar[0] = Rpar[0]*cm;
        G4double pY = Rpar[1] = Rpar[1]*cm;
        G4double pZ = Rpar[2] = Rpar[2]*cm;
        solid = new G4Box(vname, pX, pY, pZ);
    
        for (G4int i=0; i<3; i++) {
            rangehi[i] = Rpar[i];
            rangelo[i] = -Rpar[i];
        }
    }
    if ( shape == "TRD1" ) {
        G4double pdx1 = Rpar[0] = Rpar[0]*cm;
        G4double pdx2 = Rpar[1] = Rpar[1]*cm;
        G4double pdy1 = Rpar[2] = Rpar[2]*cm;
        G4double pdy2 = pdy1;
        G4double pdz  = Rpar[3] = Rpar[3]*cm;
        solid = new G4Trd(vname, pdx1, pdx2, pdy1, pdy2, pdz);
    
        rangehi[0] = max(pdx1, pdx2);
        rangelo[0] = -rangehi[0];
        rangehi[1] = pdy1;
        rangelo[1] = -rangehi[1];
        rangehi[2] = pdz;
        rangelo[2] = -rangehi[2];
    }
    if ( shape == "TRD2" ) {
        G4double pdx1 = Rpar[0] = Rpar[0]*cm;
        G4double pdx2 = Rpar[1] = Rpar[1]*cm;
        G4double pdy1 = Rpar[2] = Rpar[2]*cm;
        G4double pdy2 = Rpar[3] = Rpar[3]*cm;
        G4double pdz  = Rpar[4] = Rpar[4]*cm;
        solid = new G4Trd(vname, pdx1, pdx2, pdy1, pdy2, pdz);
    
        rangehi[0] = max(pdx1, pdx2);
        rangelo[0] = -rangehi[0];
        rangehi[1] = max(pdy1, pdy2);
        rangelo[1] = -rangehi[1];
        rangehi[2] = pdz;
        rangelo[2] = -rangehi[2];
    }
    if ( shape == "TRAP" ) {
        G4double pDz    = Rpar[0] = Rpar[0]*cm;
        G4double pTheta = Rpar[1] = Rpar[1]*deg;
        G4double pPhi   = Rpar[2] = Rpar[2]*deg;
        G4double pDy1   = Rpar[3] = Rpar[3]*cm;
        G4double pDx1   = Rpar[4] = Rpar[4]*cm;
        G4double pDx2   = Rpar[5] = Rpar[5]*cm;
        G4double pAlp1  = Rpar[6] = Rpar[6]*deg;
        G4double pDy2   = Rpar[7] = Rpar[7]*cm;
        G4double pDx3   = Rpar[8] = Rpar[8]*cm;
        G4double pDx4   = Rpar[9] = Rpar[9]*cm;
        G4double pAlp2  = Rpar[10]= Rpar[10]*deg;
        solid = new G4Trap(vname, pDz, pTheta, pPhi, pDy1, pDx1, pDx2, pAlp1,
                           pDy2, pDx3, pDx4, pAlp2);
    
            // only legal division is along z
        rangehi[0] = 0.;
        rangelo[0] = -rangehi[0];
        rangehi[1] = 0.;
        rangelo[1] = -rangehi[1];
        rangehi[2] = pDz;
        rangelo[2] = -rangehi[2];
    }
    if ( shape == "TUBE" ) {
        G4double pRMin = Rpar[0] = Rpar[0]*cm;
        G4double pRMax = Rpar[1] = Rpar[1]*cm;
        G4double pDz   = Rpar[2] = Rpar[2]*cm;
        G4double pSPhi = 0.*deg;
        G4double pDPhi = 360.*deg;
        solid = new G4Tubs(vname, pRMin, pRMax, pDz, pSPhi, pDPhi);
        rangehi[0] = pRMax;
        rangelo[0] = pRMin;
        rangehi[1] = pSPhi + pDPhi;
        rangelo[1] = pSPhi;
        rangehi[2] = pDz;
        rangelo[2] = -rangehi[2];
        axis = kRho;
    }
    if ( shape == "TUBS" ) {
        G4double pRMin = Rpar[0] = Rpar[0]*cm;
        G4double pRMax = Rpar[1] = Rpar[1]*cm;
        G4double pDz   = Rpar[2] = Rpar[2]*cm;
        G4double pSPhi = Rpar[3] = Rpar[3]*deg;
        Rpar[4] = Rpar[4]*deg;
        G4double pDPhi = Rpar[4] - pSPhi;
        if ( Rpar[4] <= pSPhi ) pDPhi = pDPhi + 360.*deg;

        solid = new G4Tubs(vname, pRMin, pRMax, pDz, pSPhi, pDPhi);
        rangehi[0] = pRMax;
        rangelo[0] = pRMin;
        rangehi[1] = pSPhi + pDPhi;
        rangelo[1] = pSPhi;
        rangehi[2] = pDz;
        rangelo[2] = -rangehi[2];
        axis = kRho;
    }
    if ( shape == "CONE" ) {
        G4double pDz    = Rpar[0] = Rpar[0]*cm;
        G4double pRmin1 = Rpar[1] = Rpar[1]*cm;
        G4double pRmax1 = Rpar[2] = Rpar[2]*cm;
        G4double pRmin2 = Rpar[3] = Rpar[3]*cm;
        G4double pRmax2 = Rpar[4] = Rpar[4]*cm;
        G4double pSPhi = 0.*deg;
        G4double pDPhi = 360.*deg;
        solid = new G4Cons(vname, pRmin1, pRmax1, pRmin2, pRmax2, pDz, pSPhi, pDPhi);
    
        rangehi[0] = max(pRmax1,pRmax2);
        rangelo[0] = min(pRmin1,pRmin2);
        rangehi[1] = pSPhi + pDPhi;
        rangelo[1] = pSPhi;
        rangehi[2] = pDz;
        rangelo[2] = -pDz;
        axis = kRho;
    }
    if ( shape == "CONS" ) {
        G4double pDz    = Rpar[0] = Rpar[0]*cm;
        G4double pRmin1 = Rpar[1] = Rpar[1]*cm;
        G4double pRmax1 = Rpar[2] = Rpar[2]*cm;
        G4double pRmin2 = Rpar[3] = Rpar[3]*cm;
        G4double pRmax2 = Rpar[4] = Rpar[4]*cm;
        G4double pSPhi  = Rpar[5] = Rpar[5]*deg;
        Rpar[6] = Rpar[6]*deg;
        G4double pDPhi  = Rpar[6]-pSPhi;
        if ( Rpar[6] <= pSPhi ) pDPhi = pDPhi + 360.*deg;

        solid = new G4Cons(vname, pRmin1, pRmax1, pRmin2, pRmax2, pDz, pSPhi, pDPhi);
    
        rangehi[0] = max(pRmax1,pRmax2);
        rangelo[0] = min(pRmin1,pRmin2);
        rangehi[1] = pSPhi + pDPhi;
        rangelo[1] = pSPhi;
        rangehi[2] = pDz;
        rangelo[2] = -pDz;
        axis = kRho;
    }
    if ( shape == "SPHE" ) {
        
        G4double pRmin  = Rpar[0] = Rpar[0]*cm;
        G4double pRmax  = Rpar[1] = Rpar[1]*cm;
        G4double pThe1  = Rpar[2] = Rpar[2]*deg;
        G4double pThe2  = Rpar[3] = Rpar[3]*deg;
        G4double pDThe  = pThe2 - pThe1;
        G4double pPhi1  = Rpar[4] = Rpar[4]*deg;
        G4double pPhi2  = Rpar[5] = Rpar[5]*deg;
        G4double pDPhi  = pPhi2 - pPhi1;
        
        rangehi[0] = pRmax;
        rangelo[0] = pRmin;
        rangehi[1] = max(pThe1, pThe2);
        rangelo[1] = min(pThe1, pThe2);
        rangehi[2] = max(pPhi1, pPhi2);
        rangehi[2] = min(pPhi1, pPhi2);
        
        solid = new G4Sphere(vname, pRmin, pRmax, pThe1, pDThe, pPhi1, pDPhi);
            
        axis = kRadial3D;
    }
    if ( shape == "PARA" ) {
        G4double pDx = Rpar[0]*cm;
        G4double pDy = Rpar[1]*cm;
        G4double pDz = Rpar[2]*cm;
        G4double pAlph = Rpar[3]*deg;
        G4double pThet = Rpar[4]*deg;
        G4double pPhi  = Rpar[5]*deg;

        solid = new G4Para(vname, pDx, pDy, pDz, pAlph, pThet, pPhi);

            // ranges given below are only correct for boxes.
        for (G4int i=0; i<2; i++) {
            rangehi[i] = Rpar[i];
            rangelo[i] = -rangehi[i];
        }
    }
    if ( shape == "PGON" ) {
        G4int i;
        G4int npdv = int(Rpar[2]);
        G4int nz = int(Rpar[3]);

        G4double pPhi1 = Rpar[0]*deg;
        G4double dPhi  = Rpar[1]*deg;
    
        G4double *DzArray = new G4double[nz];
        G4double *Rmax    = new G4double[nz];
        G4double *Rmin    = new G4double[nz];

        rangehi[0] = -kInfinity  ;
        rangelo[0] =  kInfinity ;
        rangehi[2] = -kInfinity ;
        rangelo[2] =  kInfinity ;
        G4double RMIN=1.e-4*cm;
        

        for(i=0; i<nz; i++) 
        {
            int i4=3*i+4;
            int i5=i4+1;
            int i6=i4+2;
            
            DzArray[i] = Rpar[i4]*cm;
            G4String smin="G4makevol::PGON::";
            smin+=vname;
            G4String smax=smin;
            smin+="::Rmin";
            smax+="::Rmax";
            Rmin[i] = G3Bound(smin, RMIN, kInfinity, Rpar[i5]*cm);
            Rmax[i] = G3Bound(smax, Rmin[i]+RMIN, kInfinity,
                              Rpar[i6]*cm);
            rangelo[0] = min(rangelo[0], Rmin[i]);
            rangehi[0] = max(rangehi[0], Rmax[i]);
            rangelo[2] = min(rangelo[2], DzArray[i]);
            rangehi[2] = max(rangehi[2], DzArray[i]);
        }
        for (i=0;i<nz;i++){
            assert(Rmin[i]>=0 && Rmax[i]>=Rmin[i]);
        }
        solid = new G4BREPSolidPolyhedra(vname, pPhi1, dPhi, npdv, nz,
                                         DzArray[0], DzArray, Rmin, Rmax);

        rangehi[1] = pPhi1 + dPhi;
        rangelo[1] = pPhi1;
        axis = kRho;

        delete [] DzArray;
        delete [] Rmin;
        delete [] Rmax;

    }
    if ( shape == "PCON" ) {

        G4int i;
        G4double pPhi1 = Rpar[0] = Rpar[0]*deg;
        G4double dPhi  = Rpar[1] = Rpar[1]*deg;    
        G4int nz = int(Rpar[2]);
    
        G4double *DzArray = new G4double[nz];
        G4double *Rmax    = new G4double[nz];
        G4double *Rmin    = new G4double[nz];

        rangehi[0] = -kInfinity  ;
        rangelo[0] =  kInfinity ;
        rangehi[2] = -kInfinity ;
        rangelo[2] =  kInfinity ;

        G4double RMIN=1.e-4*cm;
        
        for(i=0; i<nz; i++){
            int i4=3*i+3;
            int i5=i4+1;
            int i6=i4+2;
            
            DzArray[i] = Rpar[i4]*cm;
            G4String smin="G4makevol::PCON::";
            smin+=vname;
            G4String smax=smin;
            smin+="::Rmin";
            smax+="::Rmax";
//            Rmin[i] = G3Bound(smin, RMIN, kInfinity, Rpar[i5]*cm);
            Rmin[i] = G3Bound(smin, 0., kInfinity, Rpar[i5]*cm);
            Rmax[i] = G3Bound(smax, Rmin[i]+RMIN, kInfinity, Rpar[i6]*cm);
            rangelo[0] = min(rangelo[0], Rmin[i]);
            rangehi[0] = max(rangehi[0], Rmax[i]);
            rangelo[2] = min(rangelo[2], DzArray[i]);
            rangehi[2] = max(rangehi[2], DzArray[i]);
        }
        for (i=0;i<nz;i++){
            assert(Rmin[i]>=0 && Rmax[i]>=Rmin[i]);
        }
        solid = new G4BREPSolidPCone(vname, pPhi1, dPhi, nz,
                                     DzArray[0], DzArray, Rmin, Rmax);
        rangehi[1] = pPhi1 + dPhi;
        rangelo[1] = pPhi1;
        axis = kRho;

        delete [] DzArray;
        delete [] Rmin;
        delete [] Rmax;
    }
    if ( shape == "ELTU" ) {
            // $$$ not implemented.
            // $$$        solid = new G4Eltu(vname, Rpar);
    }
    if ( shape == "HYPE" ) {
            // $$$ not implemented.
            // $$$        solid = new G4Hype(vname, Rpar, kDegrees);
    }
    if ( shape == "GTRA" ) {
            // $$$ not implemented.
            // $$$        solid = new G4Gtra(vname, Rpar, kDegrees);
    }
    if ( shape == "CTUB" ) {
            // $$$ not implemented.
            // $$$        solid = new G4Ctub(vname, Rpar, kDegrees);
    }
  
        // get the material corresponding to the tracking medium
    G4Material *material;
        // get the magnetic field and user limits
    G4MagneticField *field;
    G4UserLimits *limits;
    G4VSensitiveDetector *sensitive = NULL;
    if (nmed>0) {
        material = G3Med.GetMat(nmed);
        field = G3Med.GetMag(nmed);
        limits = G3Med.GetLim(nmed);
    } else {
        material = NULL;
        field = NULL;
        limits = NULL;
    }

    G4LogicalVolume* lvol;
    
        // Create the logical volume
    if ( solid == NULL) {
        G4cout << "For volume " << vname << ", shape " << shape
             << " not supported " << endl;
        lvol = NULL;
    } else {

            // the GetLVx routine queries the G4LogicalVolumeStore to Retrieve
            // logical volume pointers by name
        
        lvol = G3Vol.GetLVx(vname);

        if (!lvol)
        {
	  G4FieldManager* FieldMgr = new G4FieldManager(field);
	  lvol = new G4LogicalVolume(solid, material, vname,
				     FieldMgr, sensitive, limits);
            
                // Store the logical volume name/pointer association
            G3Vol.PutLV(&vname,lvol,nmed,rangehi,rangelo,axis,shape,Rpar,npar,
                        solid);
        }
        else{
            G4cout << "G3makevol: LV " << vname << " already defined" << endl;
        }

            // check that the G3toG4 global mother is set
        
        G3Vol.SetMother(lvol);
        
            // set the visual attributes, A. Mokhtarani 2-5-97

        lvol->SetVisAttributes(0);
    }
    return lvol;
}
