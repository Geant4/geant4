// $Id: G4FlatSurface.hh,v 1.1 2003-11-12 15:43:00 link Exp $
#ifndef __G4FLATSURFACE__
#define __G4FLATSURFACE__
//*************************************************************************
//* --------------------
//* G4FlatSurface
//* --------------------
//* (Description)
//* 	Class for flat boundary surface for J4Solid.
//*     
//* (Update Record)
//*	2002/08/01  K.Hoshina	Original version.
//*************************************************************************

#include "G4VSurface.hh"

class G4TwistedTubs;

class G4FlatSurface : public G4VSurface
{
public:

   G4FlatSurface(const G4String         &name,
                 const G4RotationMatrix &rot,
                 const G4ThreeVector    &tlate,
                 const G4ThreeVector    &n,
                 const EAxis             axis1 = kRho,
                 const EAxis             axis2 = kPhi,
                       G4double          axis0min = -kInfinity,
                       G4double          axis1min = -kInfinity,
                       G4double          axis0max = kInfinity,
                       G4double          axis1max = kInfinity );
                       
   G4FlatSurface(const G4String            &name,
                       G4TwistedTubs       *solid,
                       G4int                handedness) ;
                       
   virtual ~G4FlatSurface();
   
  virtual G4ThreeVector  GetNormal(const G4ThreeVector & /* xx */ , G4bool isGlobal = FALSE);

   virtual G4int          DistanceToSurface(const G4ThreeVector &gp,
                                            const G4ThreeVector &gv,
                                                  G4ThreeVector  gxx[],
                                                  G4double       distance[],
                                                  G4int          areacode[],
                                                  G4bool         isvalid[],
                                                  EValidate      validate = kValidateWithTol);
                                                  
   virtual G4int          DistanceToSurface(const G4ThreeVector &gp,
                                                  G4ThreeVector  gxx[],
                                                  G4double       distance[],
                                                  G4int          areacode[]);
                                                  
protected:                                                  
   virtual G4int          GetAreaCode(const G4ThreeVector &xx, 
                                      G4bool withTol = TRUE) ;

private:

   virtual void           SetCorners();
   virtual void           SetBoundaries();


private:
   
};

#endif

