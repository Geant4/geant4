// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3Division.cc,v 1.2 1999-12-05 17:50:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, V.Berejnoi 13.10.99

#include "G3Division.hh"
#include "G3VolTableEntry.hh"
#include "G3toG4MakeSolid.hh"
#include "G3Pos.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

G3VolTableEntry* G4CreateVTE(G4String vname, G4String shape, G4int nmed,
                               G4double Rpar[], G4int npar);

G3Division::G3Division(G3DivType type, G3VolTableEntry* vte, 
                G3VolTableEntry* mvte, G4int nofDivisions, 
		G4int iaxis, G4int nmed, G4double c0, G4double step)
  : fType(type),
    fVTE(vte),
    fMVTE(mvte),
    fNofDivisions(nofDivisions),
    fIAxis(iaxis),
    fNmed(nmed),
    fC0(c0),
    fStep(step),
    fLowRange(0.),
    fHighRange(0.),
    fWidth(0.),
    fOffset(0.),
    fAxis(kXAxis)
{
  fVTE->SetHasNegPars(true);
}

G3Division::G3Division(G3VolTableEntry* vte, G3VolTableEntry* mvte,
                       const G3Division& division)
  : fVTE(vte),
    fMVTE(mvte)
{
  // only "input" parameters are copied from division
  fType = division.fType;
  fNofDivisions = division.fNofDivisions;
  fIAxis = division.fIAxis;
  fNmed = division.fNmed;
  fC0 = division.fC0;
  fStep = division.fStep;

  // other parameters are set as in standard constructor
  fLowRange = 0.;
  fHighRange = 0.;
  fWidth = 0.;
  fOffset = 0.;
  fAxis = kXAxis;
  fVTE->SetHasNegPars(true);
}

G3Division::~G3Division()
{}

// public methods

void G3Division::UpdateVTE()
{
  if (fVTE->HasNegPars() && !(fMVTE->HasNegPars())) {

    // set nmed from mother 
    if (fNmed == 0) fNmed = fMVTE->GetNmed();
    fVTE->SetNmed(fNmed);
   
    SetRangeAndAxis();
    
    // create envelope (if necessary)
    // and solid    
    G3VolTableEntry* envVTE = 0;
    if      (fType == kDvn)  envVTE = Dvn(); 
    else if (fType == kDvn2) envVTE = Dvn2(); 
    else if (fType == kDvt)  envVTE = Dvt(); 
    else if (fType == kDvt2) envVTE = Dvt2();
    
    if (envVTE) {
      // reset mother <-> daughter
       fMVTE->ReplaceDaughter(fVTE, envVTE);
      fVTE->ReplaceMother(fMVTE, envVTE);
      envVTE->AddDaughter(fVTE);
      envVTE->AddMother(fMVTE);

      // replace mother with envelope
      fMVTE = envVTE;
    }
  }  
}

G4VPhysicalVolume* G3Division::CreatePVReplica()
{
  G4String name = fVTE->GetName();
  G4LogicalVolume* lv =  fVTE->GetLV();
  G4LogicalVolume* mlv = fMVTE->GetLV();
  
  G4PVReplica* pvol 
    = new G4PVReplica(name, lv, mlv, fAxis, fNofDivisions, fWidth, fOffset);

  return pvol;
}

// private methods

void G3Division::Exception(G4String where, G4String what) {
  G4Exception("G3Division::" + where + " for " + what + " is not implemented");
}  

void G3Division::SetRangeAndAxis()
// set fHighRange, fLowRange, fAxis
{
    G4String shape = fMVTE->GetShape();
    G4double *Rpar = fMVTE->GetRpar();
    
    switch (fIAxis) {
      case 1: fAxis = kXAxis;
              break;
      case 2: fAxis = kYAxis;
              break;
      case 3: fAxis = kZAxis;
              break;
      default: G4Exception("G3Division: wrong iaxis defenition");
    }

    if ( shape == "BOX" ) {
      fHighRange = Rpar[fIAxis-1]*cm;
      fLowRange = -fHighRange;
    }
    else if ( shape == "TRD1" ) {
      if (fIAxis == 1){
        fHighRange = max(Rpar[0]*cm, Rpar[1]*cm);
      }
      else if( fIAxis == 2) {
       fHighRange = Rpar[2]*cm;
      }
      else if( fIAxis == 3) {
       fHighRange = Rpar[3]*cm;
      }
      fLowRange = - fHighRange;
    }
    else if ( shape == "TRD2" ) {
      if (fIAxis == 1){
        fHighRange = max(Rpar[0]*cm, Rpar[1]*cm);
      }
      else if( fIAxis == 2) {
        fHighRange = max(Rpar[2]*cm, Rpar[3]*cm);
      }
      else if( fIAxis == 3) {
       fHighRange = Rpar[4]*cm;
      }
    }
    else if ( shape == "TRAP" ) {
      if ( fIAxis == 3 ) fHighRange = Rpar[0]*cm;
      else               fHighRange = 0.;
      fLowRange = -fHighRange;
    }
    else if ( shape == "TUBE" ) {
      if (fIAxis == 1){
        fHighRange = Rpar[1]*cm;
        fLowRange = Rpar[0]*cm;
        fAxis = kRho;
      }
      else if( fIAxis == 2) {
        fHighRange = 360.*deg;
        fLowRange = 0.;
        fAxis = kPhi;
      }
      else if( fIAxis == 3) {
       fHighRange = Rpar[2]*cm;
       fLowRange = -fHighRange;
      }
    }
    else if ( shape == "TUBS" ) {
      if (fIAxis == 1){
        fHighRange = Rpar[1]*cm;
        fLowRange = Rpar[0]*cm;
        fAxis = kRho;
      }
      else if( fIAxis == 2) {

       fLowRange = Rpar[3]*deg;
       fHighRange = Rpar[4]*deg - fLowRange;
       if ( Rpar[4]*deg <= fLowRange )fHighRange = fHighRange + 360.*deg;
       fHighRange = fHighRange + fLowRange;
       fAxis = kPhi;
      }
      else if( fIAxis == 3) {
       fHighRange = Rpar[2]*cm;
       fLowRange = -fHighRange;
      }
    }
    else if ( shape == "CONE" ) {
      if (fIAxis == 1){
        fHighRange = max(Rpar[2]*cm,Rpar[4]*cm);
        fLowRange = max(Rpar[1]*cm,Rpar[3]*cm);
        fAxis = kRho;
      }
      else if( fIAxis == 2) {

       fLowRange = 0.;
       fHighRange = 360.*deg;
       fAxis = kPhi;
      }
      else if( fIAxis == 3) {
       fHighRange = Rpar[0]*cm;
       fLowRange = -fHighRange;
      }
    }
    else if ( shape == "CONS" ) {
      if (fIAxis == 1){
        fHighRange = max(Rpar[2]*cm,Rpar[4]*cm);
        fLowRange = max(Rpar[1]*cm,Rpar[3]*cm);
        fAxis = kRho;
      }
      else if( fIAxis == 2) {

       fLowRange = Rpar[5]*deg;
       fHighRange = Rpar[6]*deg - fLowRange;
       if ( Rpar[6]*deg <= fLowRange )fHighRange = fHighRange + 360.*deg;
       fHighRange = fHighRange + fLowRange;
       fAxis = kPhi;
      }
      else if( fIAxis == 3) {
       fHighRange = Rpar[2]*cm;
       fLowRange = -fHighRange;
      }
    }
    else if ( shape == "SPHE" ) {
      if (fIAxis == 1){
        fHighRange = Rpar[1]*cm;
        fLowRange = Rpar[0]*cm;
        fAxis = kRho;
      }
      else if( fIAxis == 2) {
       fLowRange = min(Rpar[2]*deg,Rpar[3]*deg);
       fHighRange = max(Rpar[2]*deg,Rpar[3]*deg);
       fAxis = kPhi;
      }
      else if( fIAxis == 3) {
       fLowRange = min(Rpar[4]*deg,Rpar[5]*deg);
       fHighRange = max(Rpar[4]*deg,Rpar[5]*deg);
       fAxis = kPhi; // ?????? 
      }
    }
    else if ( shape == "PARA" ) {
      fHighRange = Rpar[fIAxis-1]*cm;
      fLowRange = -fHighRange;
    }
    else if ( shape == "PGON" ) {
        G4int i;
        G4int npdv = int(Rpar[2]);
        G4int nz = int(Rpar[3]);

        G4double pPhi1 = Rpar[0]*deg;
        G4double dPhi  = Rpar[1]*deg;
    
        G4double *DzArray = new G4double[nz];
        G4double *Rmax    = new G4double[nz];
        G4double *Rmin    = new G4double[nz];
        G4double rangehi[3], rangelo[3];
        rangehi[0] = -kInfinity  ;
        rangelo[0] =  kInfinity ;
        rangehi[2] = -kInfinity ;
        rangelo[2] =  kInfinity ;

        for(i=0; i<nz; i++) 
        {
            int i4=3*i+4;
            int i5=i4+1;
            int i6=i4+2;
            
            DzArray[i] = Rpar[i4]*cm;
            Rmin[i] = Rpar[i5]*cm;
            Rmax[i] = Rpar[i6]*cm;
            rangelo[0] = min(rangelo[0], Rmin[i]);
            rangehi[0] = max(rangehi[0], Rmax[i]);
            rangelo[2] = min(rangelo[2], DzArray[i]);
            rangehi[2] = max(rangehi[2], DzArray[i]);
        }
        for (i=0;i<nz;i++){
            assert(Rmin[i]>=0 && Rmax[i]>=Rmin[i]);
        }
        rangehi[1] = pPhi1 + dPhi;
        rangelo[1] = pPhi1;
        fHighRange = rangehi[fIAxis-1];
        fLowRange = rangelo[fIAxis-1];
        if      (fIAxis == 1)fAxis = kRho;
        else if (fIAxis == 2)fAxis = kPhi;
        else if (fIAxis == 3)fAxis = kZAxis;

        delete [] DzArray;
        delete [] Rmin;
        delete [] Rmax;

    }
    else if ( shape == "PCON" ) {

        G4int i;
        G4double pPhi1 = Rpar[0]*deg;
        G4double dPhi  = Rpar[1]*deg;    
        G4int nz = int(Rpar[2]);
    
        G4double *DzArray = new G4double[nz];
        G4double *Rmax    = new G4double[nz];
        G4double *Rmin    = new G4double[nz];
        G4double rangehi[3],rangelo[3];

        rangehi[0] = -kInfinity  ;
        rangelo[0] =  kInfinity ;
        rangehi[2] = -kInfinity ;
        rangelo[2] =  kInfinity ;
        
        for(i=0; i<nz; i++){
            int i4=3*i+3;
            int i5=i4+1;
            int i6=i4+2;
            
            DzArray[i] = Rpar[i4]*cm;
            Rmin[i] = Rpar[i5]*cm;
            Rmax[i] = Rpar[i6]*cm;
            rangelo[0] = min(rangelo[0], Rmin[i]);
            rangehi[0] = max(rangehi[0], Rmax[i]);
            rangelo[2] = min(rangelo[2], DzArray[i]);
            rangehi[2] = max(rangehi[2], DzArray[i]);
        }
        for (i=0;i<nz;i++){
            assert(Rmin[i]>=0 && Rmax[i]>=Rmin[i]);
        }
        rangehi[1] = pPhi1 + dPhi;
        rangelo[1] = pPhi1;
        fHighRange = rangehi[fIAxis-1];
        fLowRange = rangelo[fIAxis-1];
        if      (fIAxis == 1)fAxis = kRho;
        else if (fIAxis == 2)fAxis = kPhi;
        else if (fIAxis == 3)fAxis = kZAxis;


        delete [] DzArray;
        delete [] Rmin;
        delete [] Rmax;
    }
    else if ( shape == "ELTU" ||  shape == "HYPE" || shape == "GTRA" ||
         shape == "CTUB") {
       Exception("SetRangeAndAxis", shape);
    }
    else {
       Exception("SetRangeAndAxis", "Unknown shape" + shape);
    }  

    // verbose
    //G4cout << "SetRangeAndAxis: " 
    //       << fLowRange << " " << fHighRange << " " << fAxis << endl;
}

G3VolTableEntry* G3Division::CreateEnvelope(G4String shape, G4double hi, 
                               G4double lo, G4double par[], G4int npar)
// create new VTE with G3Pos corresponding to the
// envelope of divided volume
{
    // verbose
    // G4cout << "  G3Division::CreateEnvelope " << "fIAaxis= " << fIAxis
    //        << " hi= " << hi
    //        << " lo= " << lo
    //        << endl;

    G4double *Rpar = new G4double[npar+2];
    for (G4int i=0; i<npar; ++i){ Rpar[i] = par[i];}
    G4double pos[3] = {0.,0.,0.};
  
    if ( shape == "BOX" ) {
      Rpar[fIAxis-1] = (hi - lo)/2./cm;
      pos [fIAxis-1] = (hi + lo)/2.;
    }
    else if ( shape == "TRD1" ) {
      if ( fIAxis == 1 || fIAxis == 2  ) {
        Exception("CreateEnvelope","TRD1-x,y");
      }
      else if ( fIAxis == 3 ) {
	// x = x1 + (c-z1)(x2 -x1)/(z2-z1)
	G4double tn, x1, z1;
        tn = (Rpar[1] - Rpar[0])/(2.* Rpar[3]); 
        x1 = Rpar[0]; z1 = -Rpar[3];
        Rpar[0] = x1 + tn * (lo/cm - z1);
        Rpar[1] = x1 + tn * (hi/cm - z1);
        Rpar[3] = (hi - lo)/2./cm;
        pos[2]  = (hi + lo)/2.;
      }
    }
    else if ( shape == "TRD2" ) {
      if ( fIAxis == 1 || fIAxis == 2) {
        Exception("CreateEnvelope","TRD2-x,y");
      }
      else if ( fIAxis == 3 ) {
	// x = x1 + (c-z1)(x2 -x1)/(z2-z1)
	// y = y1 + (c-z1)(y2 -y1)/(z2-z1)
	G4double tn1, tn2, x1, y1, z1;
        tn1 = (Rpar[1] - Rpar[0])/(2.* Rpar[4]); 
        tn2 = (Rpar[3] - Rpar[2])/(2.* Rpar[4]); 
        x1 = Rpar[0]; y1 = Rpar[2]; z1 = -Rpar[3];
        Rpar[0] = x1 + tn1 * (lo/cm - z1);
        Rpar[1] = x1 + tn1 * (hi/cm - z1);
        Rpar[2] = y1 + tn2 * (lo/cm - z1);
        Rpar[3] = y1 + tn2 * (hi/cm - z1);
        Rpar[4] = (hi - lo)/2./cm;
        pos[2]  = (hi + lo)/2.;
      }
    }
    else if ( shape == "TRAP" ) {
      Exception("CreateEnvelope","TRAP-x,y,z");
    }
    else if ( shape == "TUBE" ) {
      if ( fIAxis == 1 ) {
        Rpar[0] = lo/cm;
        Rpar[1] = hi/cm;
      }
      else if ( fIAxis == 2 ) {
        Rpar[3] = lo/deg;
        Rpar[4] = hi/deg;
        npar = npar + 2;
        shape = "TUBS";
      }
      else if ( fIAxis == 3 ) {
        Rpar[2] = (hi - lo)/2./cm;
        pos [2] = (hi + lo)/2.;
      }
    }
    else if ( shape == "TUBS" ) {
      if ( fIAxis == 1 ) {
        Rpar[0] = lo/cm;
        Rpar[1] = hi/cm;
      }
      else if ( fIAxis == 2 ) {
        Rpar[3] = lo/deg;
        Rpar[4] = hi/deg;
      }
      else if ( fIAxis == 3 ) {
        Rpar[2] = (hi - lo)/2./cm;
        pos [2] = (hi + lo)/2.;
      }
    }
    else if ( shape == "CONE" ) {
      if ( fIAxis == 1) {
        Exception("CreateEnvelope","CONE-x,z");
      }
      else if ( fIAxis == 2 ) {
        Rpar[5] = lo/deg;
        Rpar[6] = hi/deg;
        npar = npar + 2;
        shape = "CONS";
      }
      else if ( fIAxis == 3 ) {
        G4double tn1, tn2, rmin, rmax, z1;
        tn1 = (Rpar[3] - Rpar[1])/(2.* Rpar[0]); 
        tn2 = (Rpar[4] - Rpar[2])/(2.* Rpar[0]); 
        rmin = Rpar[1]; rmax = Rpar[2]; z1 = -Rpar[0];
        Rpar[1] = rmin + tn1 * (lo/cm - z1);
        Rpar[3] = rmin + tn1 * (hi/cm - z1);
        Rpar[2] = rmax + tn2 * (lo/cm - z1);
        Rpar[4] = rmax + tn2 * (hi/cm - z1);
        Rpar[0] = (hi - lo)/2./cm;
        pos[2]  = (hi + lo)/2.;
      }
    }
    else if ( shape == "CONS" ) {
      if ( fIAxis == 1 ) {
        Exception("CreateEnvelope","CONS-x");
      }
      else if ( fIAxis == 2 ) {
        Rpar[5] = lo/deg;
        Rpar[6] = hi/deg;
      }
      else if ( fIAxis == 3 ) {
        G4double tn1, tn2, rmin, rmax, z1;
        tn1 = (Rpar[3] - Rpar[1])/(2.* Rpar[0]); 
        tn2 = (Rpar[4] - Rpar[2])/(2.* Rpar[0]); 
        rmin = Rpar[1]; rmax = Rpar[2]; z1 = -Rpar[0];
        Rpar[1] = rmin + tn1 * (lo/cm - z1);
        Rpar[3] = rmin + tn1 * (hi/cm - z1);
        Rpar[2] = rmax + tn2 * (lo/cm - z1);
        Rpar[4] = rmax + tn2 * (hi/cm - z1);
        Rpar[0] = (hi - lo)/2./cm;
        pos[2]  = (hi + lo)/2.;
      }
    }
    else if ( shape == "SPHE" ) {
      Exception("CreateEnvelope","SPHE-x,y,z");                
    }
    else if ( shape == "PARA" ) {
      Exception("CreateEnvelope","PARA-x,y,z");
    }
    else if ( shape == "PGON" ) {
      if ( fIAxis == 2) {
	Rpar[0] = lo/deg;
	Rpar[1] = hi/deg;
	// rotm = ???
      }
      else {
        Exception("CreateEnvelope","PGON-x,z");
      }
    }
    else if ( shape == "PCON" ) {
      if ( fIAxis == 2) {
	Rpar[0] = lo/deg;
	Rpar[1] = hi/deg;
	// rotm = ???
      }
      else {
        Exception("CreateEnvelope","PCON-x,z");
      }
    }
    else {
       Exception("CreateEnvelope", "Unknown shape" + shape);
    }  

    // create new VTE corresponding to envelope
    G4String envName = fVTE->GetName() + "_ENV"; 
    G3VolTableEntry* envVTE 
      = G4CreateVTE(envName, shape, fNmed, Rpar, npar);

    // create a G3Pos object and add it to envVTE
    G4String motherName = fMVTE->GetMasterClone()->GetName();
    G4ThreeVector* offset = new G4ThreeVector(pos[0],pos[1],pos[2]);    
    G4String only = "ONLY";
    G3Pos* aG3Pos = new G3Pos(motherName, 1, offset, 0, only);              
    envVTE->AddG3Pos(aG3Pos);

    delete [] Rpar; 

    return envVTE;
}

void G3Division::CreateSolid(G4String shape, G4double par[], G4int npar)
// create the solid corresponding to divided volume
// and set the fOffset for replica
{
    G4double *Rpar = new G4double[npar+2];
    for (G4int i=0; i<npar; ++i){ Rpar[i] = par[i];}

    // verbose
    // G4cout << "G3Division::CreateSolid volume before: " 
    //        << fVTE->GetName() << " " << shape << endl;    
    // G4cout << " npar,Rpar: " << npar;
    // for (G4int ii = 0; ii < npar; ++ii) G4cout << " " << Rpar[ii];
    // G4cout << endl;
  
    if ( shape == "BOX" ) {
      if      ( fIAxis == 1 ) Rpar[0] = fWidth/2./cm;
      else if ( fIAxis == 2 ) Rpar[1] = fWidth/2./cm; 
      else if ( fIAxis == 3 ) Rpar[2] = fWidth/2./cm; 
    }
    else if ( shape == "TRD1" ) {
      if ( fIAxis == 1 || fIAxis == 2 ) {
        Exception("CreateSolid", "TRD1-x,y");
      }
      else if ( fIAxis == 3 ) {
         Rpar[3] = fWidth/2./cm; 
      }
    }
    else if ( shape == "TRD2" ) {
      if ( fIAxis == 1 || fIAxis == 2 ) {
        Exception("CreateSolid", "TRD2-x,y");
      }
      else if ( fIAxis == 3 ) {
         Rpar[4] =  fWidth/2./cm; 
      }
    }
    else if ( shape == "TRAP" ) {
      if ( fIAxis == 1 || fIAxis == 2) {
        Exception("CreateSolid", "TRAP-x,y");
      }
      else if ( fIAxis == 3 ) {
         Rpar[0] =  fWidth/2./cm; 
      }
    }
    else if ( shape == "TUBE" ) {
      if ( fIAxis == 1 ) {
         Rpar[1] = Rpar[0] + fWidth/cm;
         fOffset = Rpar[0]*cm;
      }
      else if ( fIAxis == 2 ) {
         Rpar[3] = 0.; 
         Rpar[4] = fWidth/deg; 
         shape = "TUBS";
         npar = npar + 2;
      }
      else if ( fIAxis == 3 ) {
         Rpar[2] = fWidth/2./cm; 
      }
    }
    else if ( shape == "TUBS" ) {
      if ( fIAxis == 1 ) {
        Rpar[1] = Rpar[0] + fWidth/cm;
      }
      else if ( fIAxis == 2 ) {
         fOffset = Rpar[3]*deg; 
         Rpar[3] = 0.;
         Rpar[4] =  fWidth/deg;
      }
      else if ( fIAxis == 3 ) {
         Rpar[2] = fWidth/2./cm; 
      }
    }
    else if ( shape == "CONE" ) {
      if ( fIAxis == 1 ) {
        Exception("CreateSolid", "CONE-x"); 
      }
      else if ( fIAxis == 2 ) {
         Rpar[5] = 0.;
         Rpar[6] = fWidth/deg;
         shape = "CONS";
         npar = npar + 2;
      }
      else if ( fIAxis == 3 ) {
         Rpar[0] = fWidth/2./cm; 
      }
    }
    else if ( shape == "CONS" ) {
      if ( fIAxis == 1 ) {
        Exception("CreateSolid", "CONS-x"); 
      }
      else if ( fIAxis == 2 ) {
         fOffset = Rpar[5]*deg;
         Rpar[5] = 0.;
         Rpar[6] = fWidth/deg;
      }
      else if ( fIAxis == 3 ) {
         Rpar[0] = fWidth/2./cm; 
      }
    }
    else if (shape == "SPHE" || shape == "PARA") {
      Exception("CreateSolid", shape);
    }
    else if ( shape == "PGON" ) {
      if ( fIAxis == 2 ) {
         fOffset = Rpar[0]*deg;
         Rpar[0] = 0.;
         Rpar[1] = fWidth/deg;
         Rpar[2] = 1.;
      }
      else
       Exception("CreateSolid", shape);
    }
    else if ( shape == "PCON" ) {
      if ( fIAxis == 2 ) {
         fOffset = Rpar[0]*deg;
         Rpar[0] = 0.;
         Rpar[1] = fWidth/deg;
      }
      else {
        Exception("CreateSolid", shape);
      }	
    }
    else {
       Exception("CreateSolid", "Unknown shape" + shape);
    }  

    // create solid and set it to fVTE
    G4bool hasNegPars;
    G4bool deferred;   
    G4bool okAxis[3];
    G4VSolid* solid
    = G3toG4MakeSolid(fVTE->GetName(), shape, Rpar, npar, hasNegPars, deferred, okAxis);  

    if (hasNegPars) {
       G4String name = fVTE->GetName();
       G4Exception("CreateSolid VTE " + name + " has negative parameters.");
    }   
    
    // update vte
    fVTE->SetSolid(solid);
    fVTE->SetNRpar(npar, Rpar); 
    fVTE->SetHasNegPars(hasNegPars);

    // verbose
    // G4cout << "G3Division::CreateSolid volume after: " 
    //        << fVTE->GetName() << " " << shape << endl;    
    // G4cout << " npar,Rpar: " << npar;
    // for (G4int iii = 0; iii < npar; ++iii) G4cout << " " << Rpar[iii];
    // G4cout << endl;
}


G3VolTableEntry* G3Division::Dvn()
{   
  // no envelope need to be created 

  // get parameters from mother
  G4String shape = fMVTE->GetShape(); 
  G4double* Rpar = fMVTE->GetRpar();
  G4int     npar = fMVTE->GetNpar();

  // set width for replica and create solid
  fWidth = (fHighRange - fLowRange)/fNofDivisions;
  CreateSolid(shape, Rpar, npar);

  return 0;			  
}

G3VolTableEntry* G3Division::Dvn2()
{
  // to be defined as const of this class
  G4double Rmin = 0.0001*cm;

  G4String shape = fMVTE->GetShape();
  G4double* Rpar = fMVTE->GetRpar();
  G4int     npar = fMVTE->GetNpar();

  G4double c0 = fC0;
  if (fAxis == kPhi)  c0 = c0*deg;
  else                c0 = c0*cm;
          
  // create envelope (if needed)
  G3VolTableEntry* envVTE = 0;
  if( abs(c0 - fLowRange) > Rmin) {
    envVTE = CreateEnvelope(shape, fHighRange, c0, Rpar, npar);
    Rpar = envVTE->GetRpar();
    npar = envVTE->GetNpar();
  }  

  // set width for replica and create solid
  fWidth = (fHighRange - c0)/fNofDivisions;
  CreateSolid(shape, Rpar, npar);

  return envVTE;
}

G3VolTableEntry* G3Division::Dvt()
{
  // to be defined as const of this class
  G4double Rmin = 0.0001*cm;

  // get parameters from mother
  G4String shape = fMVTE->GetShape();
  G4double* Rpar = fMVTE->GetRpar();
  G4int     npar = fMVTE->GetNpar();

  // calculate the number of divisions    
  G4int ndvmx = fNofDivisions;
  G4double step = fStep;
  
  if (fAxis == kPhi)  step = step*deg;
  else                step = step*cm;

  G4int ndiv = G4int((fHighRange - fLowRange + Rmin)/step);
  // to be added warning
  if (ndvmx > 255) ndvmx = 255;
  if (ndiv > ndvmx && ndvmx > 0 ) ndiv = ndvmx;

  // create envVTE (if needed)
  G3VolTableEntry* envVTE = 0;
  G4double delta = abs((fHighRange - fLowRange) - ndiv*step);
  if (delta > Rmin) {
    envVTE 
       = CreateEnvelope(shape, fHighRange-delta/2., fLowRange+delta/2., 
                        Rpar, npar);
    Rpar = envVTE->GetRpar();
    npar = envVTE->GetNpar();
  }

  // set width for replica and create solid
  fWidth = step;
  fNofDivisions = ndiv;
  CreateSolid(shape, Rpar, npar);

  return envVTE;
}

G3VolTableEntry* G3Division::Dvt2()
{
  // to be defined as const of this class
  G4double Rmin = 0.0001*cm;

  // get parameters from mother
  G4String shape = fMVTE->GetShape();
  G4double* Rpar = fMVTE->GetRpar();
  G4int     npar = fMVTE->GetNpar();

  // calculate the number of divisions   
  G4int ndvmx = fNofDivisions;
  G4double step = fStep;
  G4double c0 = fC0;

  if(fAxis == kPhi){
    step = step*deg;
    c0 = c0*deg;
  } 
  else {
    step = step*cm;
    c0 = c0*cm;
  }  

  G4int ndiv = G4int((fHighRange - c0 + Rmin)/step);
  // to be added warning
  if (ndvmx > 255) ndvmx = 255;
  if (ndiv > ndvmx && ndvmx > 0 ) ndiv = ndvmx;

  // create envelope (if needed)
  G3VolTableEntry* envVTE = 0;
  G4double delta = abs((fHighRange - c0) - ndiv*step);
  if (abs(c0 - fLowRange) > Rmin) {
    envVTE 
      = CreateEnvelope(shape, fHighRange-delta/2., c0+delta/2., Rpar, npar);
    Rpar = envVTE->GetRpar();
    npar = envVTE->GetNpar();
  }

  // set with for replica and create solid
  fWidth = step;
  fNofDivisions = ndiv;
  CreateSolid(shape, Rpar, npar);

  return envVTE;	 
}
