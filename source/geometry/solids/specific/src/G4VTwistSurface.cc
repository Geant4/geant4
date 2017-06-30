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
// $Id: G4VTwistSurface.cc 104105 2017-05-11 08:23:18Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4VTwistSurface.cc
//
// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------

#include <iomanip>

#include "G4VTwistSurface.hh"
#include "G4GeometryTolerance.hh"

const G4int  G4VTwistSurface::sOutside        = 0x00000000;
const G4int  G4VTwistSurface::sInside         = 0x10000000;
const G4int  G4VTwistSurface::sBoundary       = 0x20000000;
const G4int  G4VTwistSurface::sCorner         = 0x40000000;
const G4int  G4VTwistSurface::sC0Min1Min      = 0x40000101; 
const G4int  G4VTwistSurface::sC0Max1Min      = 0x40000201;
const G4int  G4VTwistSurface::sC0Max1Max      = 0x40000202; 
const G4int  G4VTwistSurface::sC0Min1Max      = 0x40000102; 
const G4int  G4VTwistSurface::sAxisMin        = 0x00000101; 
const G4int  G4VTwistSurface::sAxisMax        = 0x00000202; 
const G4int  G4VTwistSurface::sAxisX          = 0x00000404;
const G4int  G4VTwistSurface::sAxisY          = 0x00000808;
const G4int  G4VTwistSurface::sAxisZ          = 0x00000C0C;
const G4int  G4VTwistSurface::sAxisRho        = 0x00001010;
const G4int  G4VTwistSurface::sAxisPhi        = 0x00001414;

// mask
const G4int  G4VTwistSurface::sAxis0          = 0x0000FF00;
const G4int  G4VTwistSurface::sAxis1          = 0x000000FF;
const G4int  G4VTwistSurface::sSizeMask       = 0x00000303;
const G4int  G4VTwistSurface::sAxisMask       = 0x0000FCFC;
const G4int  G4VTwistSurface::sAreaMask       = 0XF0000000;

//=====================================================================
//* constructors ------------------------------------------------------

G4VTwistSurface::G4VTwistSurface(const G4String &name)
  : fIsValidNorm(false), fName(name)
{

   fAxis[0]    = kUndefined;
   fAxis[1]    = kUndefined;
   fAxisMin[0] = kInfinity;
   fAxisMin[1] = kInfinity;
   fAxisMax[0] = kInfinity;
   fAxisMax[1] = kInfinity;
   fHandedness = 1;

   for (G4int i=0; i<4; i++)
   {
      fCorners[i].set(kInfinity, kInfinity, kInfinity);
      fNeighbours[i] = 0;
   }

   fCurrentNormal.p.set(kInfinity, kInfinity, kInfinity);
   
   fAmIOnLeftSide.me.set(kInfinity, kInfinity, kInfinity);
   fAmIOnLeftSide.vec.set(kInfinity, kInfinity, kInfinity);
   kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

G4VTwistSurface::G4VTwistSurface(const G4String         &name,
                       const G4RotationMatrix &rot,
                       const G4ThreeVector    &tlate,
                             G4int             handedness,
                       const EAxis             axis0 ,
                       const EAxis             axis1 ,
                             G4double          axis0min,
                             G4double          axis1min,
                             G4double          axis0max,
                             G4double          axis1max )
   : fIsValidNorm(false), fName(name)
{
   fAxis[0]    = axis0;
   fAxis[1]    = axis1;
   fAxisMin[0] = axis0min;
   fAxisMin[1] = axis1min;
   fAxisMax[0] = axis0max;
   fAxisMax[1] = axis1max;
   fHandedness = handedness;
   fRot        = rot;
   fTrans      = tlate;

   for (G4int i=0; i<4; i++)
   {
      fCorners[i].set(kInfinity, kInfinity, kInfinity);
      fNeighbours[i] = 0;
   }

   fCurrentNormal.p.set(kInfinity, kInfinity, kInfinity);
   
   fAmIOnLeftSide.me.set(kInfinity, kInfinity, kInfinity);
   fAmIOnLeftSide.vec.set(kInfinity, kInfinity, kInfinity);
   kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

//=====================================================================
//* Fake default constructor ------------------------------------------

G4VTwistSurface::G4VTwistSurface( __void__& )
  : fHandedness(0), fIsValidNorm(false), kCarTolerance(0.),
    fName("")
{
   fAxis[0] = fAxis[1] = kXAxis;
   fAxisMin[0] = fAxisMin[1] = 0.;
   fAxisMax[0] = fAxisMax[1] = 0.;
   fNeighbours[0] = fNeighbours[1] = fNeighbours[2] = fNeighbours[3] = 0;
}

//=====================================================================
//* destructor --------------------------------------------------------

G4VTwistSurface::~G4VTwistSurface()
{
}

//=====================================================================
//* AmIOnLeftSide -----------------------------------------------------

G4int G4VTwistSurface::AmIOnLeftSide(const G4ThreeVector &me, 
                                     const G4ThreeVector &vec,
                                     G4bool        withtol) 
{
   // AmIOnLeftSide returns phi-location of "me"
   // (phi relation between me and vec projected on z=0 plane).
   // If "me" is on -ve-phi-side of "vec", it returns 1.
   // On the other hand, if "me" is on +ve-phi-side of "vec",
   // it returns -1.
   // (The return value represents z-coordinate of normal vector
   //  of me.cross(vec).)
   // If me is on boundary of vec, return 0.

   const G4double kAngTolerance
     = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

   G4RotationMatrix unitrot;
   const G4RotationMatrix rottol    = unitrot.rotateZ(0.5*kAngTolerance);
   const G4RotationMatrix invrottol = unitrot.rotateZ(-1.*kAngTolerance);

   if (fAmIOnLeftSide.me == me 
       && fAmIOnLeftSide.vec == vec
       && fAmIOnLeftSide.withTol == withtol) {
      return fAmIOnLeftSide.amIOnLeftSide;
   }
   
   fAmIOnLeftSide.me      = me;
   fAmIOnLeftSide.vec     = vec;
   fAmIOnLeftSide.withTol = withtol;
   
   G4ThreeVector met   = (G4ThreeVector(me.x(), me.y(), 0.)).unit();
   G4ThreeVector vect  = (G4ThreeVector(vec.x(), vec.y(), 0.)).unit();
   
   G4ThreeVector ivect = invrottol * vect;
   G4ThreeVector rvect = rottol * vect;

   G4double metcrossvect = met.x() * vect.y() - met.y() * vect.x();
   
   if (withtol) {
      if (met.x() * ivect.y() - met.y() * ivect.x() > 0 && 
          metcrossvect >= 0)  {
         fAmIOnLeftSide.amIOnLeftSide = 1;
      } else if (met.x() * rvect.y() - met.y() * rvect.x() < 0 &&
                 metcrossvect <= 0)  {
         fAmIOnLeftSide.amIOnLeftSide = -1;
      } else {
         fAmIOnLeftSide.amIOnLeftSide = 0;
      }
   } else {
      if (metcrossvect > 0) {    
         fAmIOnLeftSide.amIOnLeftSide = 1;
      } else if (metcrossvect < 0 ) {
         fAmIOnLeftSide.amIOnLeftSide = -1;
      } else {       
         fAmIOnLeftSide.amIOnLeftSide = 0;
      }
   }

#ifdef G4TWISTDEBUG
   G4cout << "         === G4VTwistSurface::AmIOnLeftSide() =============="
          << G4endl;
   G4cout << "             Name , returncode  : " << fName << " " 
                       << fAmIOnLeftSide.amIOnLeftSide <<  G4endl;
   G4cout << "             me, vec    : " << std::setprecision(14) << me 
                                          << " " << vec  << G4endl;
   G4cout << "             met, vect  : " << met << " " << vect  << G4endl;
   G4cout << "             ivec, rvec : " << ivect << " " << rvect << G4endl;
   G4cout << "             met x vect : " << metcrossvect << G4endl;
   G4cout << "             met x ivec : " << met.cross(ivect) << G4endl;
   G4cout << "             met x rvec : " << met.cross(rvect) << G4endl;
   G4cout << "         =============================================="
          << G4endl;
#endif

   return fAmIOnLeftSide.amIOnLeftSide;
}

//=====================================================================
//* DistanceToBoundary ------------------------------------------------

G4double G4VTwistSurface::DistanceToBoundary(G4int areacode,
                                        G4ThreeVector &xx,
                                        const G4ThreeVector &p) 
{
   // DistanceToBoundary 
   //
   // return distance to nearest boundary from arbitrary point p 
   // in local coodinate.
   // Argument areacode must be one of them:
   // sAxis0 & sAxisMin, sAxis0 & sAxisMax,
   // sAxis1 & sAxisMin, sAxis1 & sAxisMax.
   //

   G4ThreeVector d;    // direction vector of the boundary
   G4ThreeVector x0;   // reference point of the boundary
   G4double      dist = kInfinity;
   G4int         boundarytype;

   if (IsAxis0(areacode) && IsAxis1(areacode)) {
      std::ostringstream message;
      message << "Point is in the corner area." << G4endl
              << "        Point is in the corner area. This function returns"
              << G4endl
              << "        a direction vector of a boundary line." << G4endl
              << "        areacode = " << areacode;
      G4Exception("G4VTwistSurface::DistanceToBoundary()", "GeomSolids0003",
                  FatalException, message);
   } else if (IsAxis0(areacode) || IsAxis1(areacode)) {
      GetBoundaryParameters(areacode, d, x0, boundarytype);
      if (boundarytype == sAxisPhi) {
         G4double t = x0.getRho() / p.getRho();
         xx.set(t*p.x(), t*p.y(), x0.z());
         dist = (xx - p).mag();
      } else { 
         // linear boundary
         // sAxisX, sAxisY, sAxisZ, sAxisRho
         dist = DistanceToLine(p, x0, d, xx);
      }
   } else {
      std::ostringstream message;
      message << "Bad areacode of boundary." << G4endl
              << "        areacode = " << areacode;
      G4Exception("G4VTwistSurface::DistanceToBoundary()", "GeomSolids0003",
                  FatalException, message);
   }
   return dist;
}

//=====================================================================
//* DistanceToIn ------------------------------------------------------

G4double G4VTwistSurface::DistanceToIn(const G4ThreeVector &gp,
                                  const G4ThreeVector &gv,
                                        G4ThreeVector &gxxbest)
{
#ifdef G4TWISTDEBUG
   G4cout << " ~~~~~ G4VTwistSurface::DistanceToIn(p,v) - Start ~~~~~" << G4endl;
   G4cout << "      Name : " << fName << G4endl;
   G4cout << "      gp   : " << gp << G4endl;
   G4cout << "      gv   : " <<  gv << G4endl;
   G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
#endif
   
   G4ThreeVector gxx[G4VSURFACENXX];
   G4double      distance[G4VSURFACENXX]  ; 
   G4int         areacode[G4VSURFACENXX]  ;
   G4bool        isvalid[G4VSURFACENXX]   ; 
   
   for (G4int i = 0 ; i<G4VSURFACENXX ; i++ ) {
     distance[i] = kInfinity ;
     areacode[i] = sOutside ;
     isvalid[i] = false ;
   }

   G4double      bestdistance   = kInfinity;
#ifdef G4TWISTDEBUG
   G4int         besti          = -1;  
#endif
   G4ThreeVector bestgxx(kInfinity, kInfinity, kInfinity);

   G4int         nxx = DistanceToSurface(gp, gv, gxx, distance, areacode, 
                                         isvalid, kValidateWithTol);

   for (G4int i=0; i< nxx; i++) {

      // skip this intersection if:
      //   - invalid intersection
      //   - particle goes outword the surface

      if (!isvalid[i]) {
         // xx[i] is sOutside or distance[i] < 0
         continue;      
      }

      G4ThreeVector normal = GetNormal(gxx[i], true);

      if ((normal * gv) >= 0) {

#ifdef G4TWISTDEBUG
         G4cout << "   G4VTwistSurface::DistanceToIn(p,v): "
                << "particle goes outword the surface." << G4endl;
#endif 
         continue; 
      }
      
      //
      // accept this intersection if the intersection is inside.
      //

      if (IsInside(areacode[i])) {
         if (distance[i] < bestdistance) {
            bestdistance = distance[i];
            bestgxx = gxx[i];
#ifdef G4TWISTDEBUG
            besti   = i;
            G4cout << "   G4VTwistSurface::DistanceToIn(p,v): "
                   << " areacode sInside name, distance = "
                   << fName <<  " "<< bestdistance << G4endl;
#endif 
         }

      //
      // else, the intersection is on boundary or corner.
      //

      } else {

         G4VTwistSurface *neighbours[2];
         G4bool      isaccepted[2] = {false, false};
         G4int       nneighbours   = GetNeighbours(areacode[i], neighbours);
            
         for (G4int j=0; j< nneighbours; j++) {
            // if on corner, nneighbours = 2.
            // if on boundary, nneighbours = 1.

            G4ThreeVector tmpgxx[G4VSURFACENXX];
            G4double      tmpdist[G4VSURFACENXX] ;
            G4int         tmpareacode[G4VSURFACENXX] ;
            G4bool        tmpisvalid[G4VSURFACENXX] ;

            for (G4int l = 0 ; l<G4VSURFACENXX ; l++ ) {
              tmpdist[l] = kInfinity ;
              tmpareacode[l] = sOutside ;
              tmpisvalid[l] = false ;
            }

            G4int tmpnxx = neighbours[j]->DistanceToSurface(
                                          gp, gv, tmpgxx, tmpdist,
                                          tmpareacode, tmpisvalid,
                                          kValidateWithTol);
            G4ThreeVector neighbournormal;

            for (G4int k=0; k< tmpnxx; k++) {

               //  
               // if tmpxx[k] is valid && sInside, the final winner must
               // be neighbour surface. return kInfinity. 
               // else , choose tmpxx on same boundary of xx, then check normal 
               //  

               if (IsInside(tmpareacode[k])) {

#ifdef G4TWISTDEBUG
                  G4cout << "   G4VTwistSurface:DistanceToIn(p,v): "
                         << " intersection "<< tmpgxx[k] << G4endl
                         << "   is inside of neighbour surface of " << fName 
                         << " . returning kInfinity." << G4endl;
                  G4cout << "~~~~~ G4VTwistSurface::DistanceToIn(p,v) - return ~~~~"
                         << G4endl;
                  G4cout << "      No intersections " << G4endl; 
                  G4cout << "      Name : " << fName << G4endl; 
                  G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                         << G4endl;
#endif 
                  if (tmpisvalid[k])  return kInfinity;
                  continue;

               //  
               // if tmpxx[k] is valid && sInside, the final winner must
               // be neighbour surface. return .  
               //

               } else if (IsSameBoundary(this,areacode[i],
                                         neighbours[j], tmpareacode[k])) { 
                  // tmpxx[k] is same boundary (or corner) of xx.
                 
                  neighbournormal = neighbours[j]->GetNormal(tmpgxx[k], true);
                  if (neighbournormal * gv < 0) isaccepted[j] = true;
               }
            } 

            // if nneighbours = 1, chabge isaccepted[1] before
            // exiting neighboursurface loop.  

            if (nneighbours == 1) isaccepted[1] = true;

         } // neighboursurface loop end

         // now, we can accept xx intersection

         if (isaccepted[0] == true && isaccepted[1] == true) {
            if (distance[i] < bestdistance) {
               bestdistance = distance[i];
               gxxbest = gxx[i];
#ifdef G4TWISTDEBUG
               besti   = i;
               G4cout << "   G4VTwistSurface::DistanceToIn(p,v): "
                      << " areacode sBoundary & sBoundary distance = "
                      << fName  << " " << distance[i] << G4endl;
#endif 
            }
         }

      } // else end
   } // intersection loop end

   gxxbest = bestgxx;

#ifdef G4TWISTDEBUG
   if (besti < 0) {
      G4cout << "~~~~~ G4VTwistSurface::DistanceToIn(p,v) - return ~~~~" << G4endl;
      G4cout << "      No intersections " << G4endl; 
      G4cout << "      Name : " << fName << G4endl; 
      G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
   } else {
      G4cout << "~~~~~ G4VTwistSurface::DistanceToIn(p,v) : return ~~~~" << G4endl;
      G4cout << "      Name, i  : " << fName << " , " << besti << G4endl; 
      G4cout << "      gxx[i]   : " << gxxbest << G4endl; 
      G4cout << "      bestdist : " << bestdistance << G4endl;
      G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
   } 

#endif

   return bestdistance;
}

//=====================================================================
//* DistanceToOut(p, v) -----------------------------------------------

G4double G4VTwistSurface::DistanceToOut(const G4ThreeVector &gp,
                                   const G4ThreeVector &gv,
                                         G4ThreeVector &gxxbest)
{
#ifdef G4TWISTDEBUG
   G4cout << "~~~~~ G4VTwistSurface::DistanceToOut(p,v) - Start ~~~~" << G4endl;
   G4cout << "      Name : " << fName << G4endl;
   G4cout << "      gp   : " << gp << G4endl;
   G4cout << "      gv   : " <<  gv << G4endl;
   G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
#endif

   G4ThreeVector gxx[G4VSURFACENXX];
   G4double      distance[G4VSURFACENXX]  ; 
   G4int         areacode[G4VSURFACENXX]  ;
   G4bool        isvalid[G4VSURFACENXX]   ; 
   G4int         i;
   
   for ( i = 0 ; i<G4VSURFACENXX ; i++ )
   {
     distance[i] = kInfinity ;
     areacode[i] = sOutside ;
     isvalid[i] = false ;
   }

   G4int         nxx;
   G4double      bestdistance   = kInfinity;

   nxx = DistanceToSurface(gp, gv, gxx, distance, areacode,
                           isvalid, kValidateWithTol);

   for (i=0; i<nxx; i++) {
      if (!(isvalid[i])) {
         continue;
      }

      G4ThreeVector normal = GetNormal(gxx[i], true);
      if (normal * gv <= 0) {
         // particle goes toword inside of solid, return kInfinity
#ifdef G4TWISTDEBUG
          G4cout << "   G4VTwistSurface::DistanceToOut(p,v): normal*gv < 0, normal " 
                 << fName << " " << normal 
                 << G4endl;
#endif 
      } else {
         // gxx[i] is accepted.
         if (distance[i] < bestdistance) {
            bestdistance = distance[i];
            gxxbest = gxx[i];
         }
      } 
   }

#ifdef G4TWISTDEBUG
   if (besti < 0) {
      G4cout << "~~~~~ G4VTwistSurface::DistanceToOut(p,v) - return ~~~" << G4endl;
      G4cout << "      No intersections   " << G4endl; 
      G4cout << "      Name     : " << fName << G4endl; 
      G4cout << "      bestdist : " << bestdistance << G4endl;
      G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
   } else {
      G4cout << "~~~~~ G4VTwistSurface::DistanceToOut(p,v) : return ~~~" << G4endl;
      G4cout << "      Name, i  : " << fName << " , " << i << G4endl; 
      G4cout << "      gxx[i]   : " << gxxbest << G4endl; 
      G4cout << "      bestdist : " << bestdistance << G4endl;
      G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
   } 
#endif

   return bestdistance;
}

//=====================================================================
//* DistanceTo(p) -----------------------------------------------------

G4double G4VTwistSurface::DistanceTo(const G4ThreeVector &gp,
                                      G4ThreeVector &gxxbest)
{
#ifdef G4TWISTDEBUG
   G4cout << "~~~~~ G4VTwistSurface::DistanceTo(p) - Start ~~~~~~~~~" << G4endl;
   G4cout << "      Name : " << fName << G4endl;
   G4cout << "      gp   : " << gp << G4endl;
   G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
#endif


   G4ThreeVector gxx[G4VSURFACENXX];
   G4double      distance[G4VSURFACENXX]  ; 
   G4int         areacode[G4VSURFACENXX]  ;
   
   for (G4int i = 0 ; i<G4VSURFACENXX ; i++ ) {
     distance[i] = kInfinity ;
     areacode[i] = sOutside ;
   }


   DistanceToSurface(gp, gxx, distance, areacode);
   gxxbest = gxx[0];

#ifdef G4TWISTDEBUG
   G4cout << "~~~~~ G4VTwistSurface::DistanceTo(p) - return ~~~~~~~~" << G4endl;
   G4cout << "      Name     : " << fName << G4endl; 
   G4cout << "      gxx      : " << gxxbest << G4endl; 
   G4cout << "      bestdist : " << distance[0] << G4endl;
   G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
#endif

   return distance[0];
}

//=====================================================================
//* IsSameBoundary ----------------------------------------------------

G4bool G4VTwistSurface::IsSameBoundary(G4VTwistSurface *surface1, G4int areacode1,
                                  G4VTwistSurface *surface2, G4int areacode2 ) const
{
   //
   // IsSameBoundary
   //
   // checking tool whether two boundaries on different surfaces are same or not.
   //

   G4bool testbitmode = true;
   G4bool iscorner[2] = {IsCorner(areacode1, testbitmode), 
                         IsCorner(areacode2, testbitmode)};

   if (iscorner[0] && iscorner[1]) {
      // on corner 
      G4ThreeVector corner1 = 
           surface1->ComputeGlobalPoint(surface1->GetCorner(areacode1));
      G4ThreeVector corner2 = 
           surface2->ComputeGlobalPoint(surface2->GetCorner(areacode2));

      if ((corner1 - corner2).mag() < kCarTolerance) {
         return true;
      } else {
         return false;
      }
    
   } else if ((IsBoundary(areacode1, testbitmode) && (!iscorner[0])) &&
              (IsBoundary(areacode2, testbitmode) && (!iscorner[1]))) {
      // on boundary  
      G4ThreeVector d1, d2, ld1, ld2;
      G4ThreeVector x01, x02, lx01, lx02;
      G4int         type1, type2;
      surface1->GetBoundaryParameters(areacode1, ld1, lx01, type1);
      surface2->GetBoundaryParameters(areacode2, ld2, lx02, type2);

      x01 = surface1->ComputeGlobalPoint(lx01);
      x02 = surface2->ComputeGlobalPoint(lx02);
      d1  = surface1->ComputeGlobalDirection(ld1);
      d2  = surface2->ComputeGlobalDirection(ld2);

      if ((x01 - x02).mag() < kCarTolerance &&
          (d1 - d2).mag() < kCarTolerance) {
        return true;
      } else {
        return false;
      }

   } else {
      return false;
   }
}

//=====================================================================
//* GetBoundaryParameters ---------------------------------------------

void G4VTwistSurface::GetBoundaryParameters(const G4int         &areacode,
                                             G4ThreeVector &d,
                                             G4ThreeVector &x0,
                                             G4int         &boundarytype) const
{
   // areacode must be one of them:
   // sAxis0 & sAxisMin, sAxis0 & sAxisMax,
   // sAxis1 & sAxisMin, sAxis1 & sAxisMax.
   
   G4int i;
   for (i=0; i<4; i++) {
      if (fBoundaries[i].GetBoundaryParameters(areacode, d, x0,
                                               boundarytype)) {
         return;
      }
   }

   std::ostringstream message;
   message << "Not registered boundary." << G4endl
           << "        Boundary at areacode " << std::hex << areacode
           << std::dec << G4endl
           << "        is not registered.";
   G4Exception("G4VTwistSurface::GetBoundaryParameters()", "GeomSolids0002",
               FatalException, message);
}

//=====================================================================
//* GetBoundaryAtPZ ---------------------------------------------------

G4ThreeVector G4VTwistSurface::GetBoundaryAtPZ(G4int areacode,
                                          const G4ThreeVector &p) const
{
   // areacode must be one of them:
   // sAxis0 & sAxisMin, sAxis0 & sAxisMax,
   // sAxis1 & sAxisMin, sAxis1 & sAxisMax.

   if (areacode & sAxis0 && areacode & sAxis1) {
     std::ostringstream message;
     message << "Point is in the corner area." << G4endl
             << "        This function returns "
             << "a direction vector of a boundary line." << G4endl
             << "        areacode = " << areacode;
     G4Exception("G4VTwistSurface::GetBoundaryAtPZ()", "GeomSolids0003",
                 FatalException, message);
   }

   G4ThreeVector d;
   G4ThreeVector x0;
   G4int         boundarytype;
   G4bool        found = false;
   
   for (G4int i=0; i<4; i++) {
      if (fBoundaries[i].GetBoundaryParameters(areacode, d, x0, 
                                                boundarytype)){
         found = true;
         continue;
      }
   }

   if (!found) {
     std::ostringstream message;
     message << "Not registered boundary." << G4endl
             << "        Boundary at areacode " << areacode << G4endl
             << "        is not registered.";
     G4Exception("G4VTwistSurface::GetBoundaryAtPZ()", "GeomSolids0002",
                 FatalException, message);
   }

   if (((boundarytype & sAxisPhi) == sAxisPhi) ||
       ((boundarytype & sAxisRho) == sAxisRho)) {
     std::ostringstream message;
     message << "Not a z-depended line boundary." << G4endl
             << "        Boundary at areacode " << areacode << G4endl
             << "        is not a z-depended line.";
     G4Exception("G4VTwistSurface::GetBoundaryAtPZ()", "GeomSolids0002",
                 FatalException, message);
   }
   return ((p.z() - x0.z()) / d.z()) * d + x0;
}

//=====================================================================
//* SetCorner ---------------------------------------------------------

void G4VTwistSurface::SetCorner(G4int areacode, G4double x, G4double y, G4double z)
{
   if ((areacode & sCorner) != sCorner){
     std::ostringstream message;
     message << "Area code must represents corner." << G4endl
             << "        areacode " << areacode;
     G4Exception("G4VTwistSurface::SetCorner()", "GeomSolids0002",
                 FatalException, message);
   }

   if ((areacode & sC0Min1Min) == sC0Min1Min) {
      fCorners[0].set(x, y, z);
   } else if ((areacode & sC0Max1Min) == sC0Max1Min) {
      fCorners[1].set(x, y, z);
   } else if ((areacode & sC0Max1Max) == sC0Max1Max) {
      fCorners[2].set(x, y, z);
   } else if ((areacode & sC0Min1Max) == sC0Min1Max) {
      fCorners[3].set(x, y, z);
   }
}

//=====================================================================
//* SetBoundaryAxis ---------------------------------------------------

void G4VTwistSurface::GetBoundaryAxis(G4int areacode, EAxis axis[]) const
{
   if ((areacode & sBoundary) != sBoundary) {
     G4Exception("G4VTwistSurface::GetBoundaryAxis()", "GeomSolids0003",
                 FatalException, "Not located on a boundary!");
   }
   G4int i;
   for (i=0; i<2; i++) {

      G4int whichaxis = 0 ;
      if (i == 0) {
         whichaxis = sAxis0;
      } else if (i == 1) {
         whichaxis = sAxis1;
      }
      
      // extracted axiscode of whichaxis
      G4int axiscode = whichaxis & sAxisMask & areacode ; 
      if (axiscode) {
         if (axiscode == (whichaxis & sAxisX)) {
            axis[i] = kXAxis;
         } else if (axiscode == (whichaxis & sAxisY)) {
            axis[i] = kYAxis;
         } else if (axiscode == (whichaxis & sAxisZ)) {
            axis[i] = kZAxis;
         } else if (axiscode == (whichaxis & sAxisRho)) {
            axis[i] = kRho;
         } else if (axiscode == (whichaxis & sAxisPhi)) {
            axis[i] = kPhi;
         } else {
           std::ostringstream message;
           message << "Not supported areacode." << G4endl
                   << "        areacode " << areacode;
           G4Exception("G4VTwistSurface::GetBoundaryAxis()", "GeomSolids0001",
                       FatalException, message);
         }
      }
   }
}

//=====================================================================
//* SetBoundaryLimit --------------------------------------------------

void G4VTwistSurface::GetBoundaryLimit(G4int areacode, G4double limit[]) const
{
   if (areacode & sCorner) {
      if (areacode & sC0Min1Min) {
         limit[0] = fAxisMin[0];
         limit[1] = fAxisMin[1];
      } else if (areacode & sC0Max1Min) {
         limit[0] = fAxisMax[0];
         limit[1] = fAxisMin[1];
      } else if (areacode & sC0Max1Max) {
         limit[0] = fAxisMax[0];
         limit[1] = fAxisMax[1];
      } else if (areacode & sC0Min1Max) {
         limit[0] = fAxisMin[0];
         limit[1] = fAxisMax[1];
      }
   } else if (areacode & sBoundary) {
      if (areacode & (sAxis0 | sAxisMin)) {
         limit[0] = fAxisMin[0];
      } else if (areacode & (sAxis1 | sAxisMin)) {
         limit[0] = fAxisMin[1];
      } else if (areacode & (sAxis0 | sAxisMax)) {
         limit[0] = fAxisMax[0];
      } else if (areacode & (sAxis1 | sAxisMax)) {
         limit[0] = fAxisMax[1];
      }
   } else {
     std::ostringstream message;
     message << "Not located on a boundary!" << G4endl
             << "          areacode " << areacode;
     G4Exception("G4VTwistSurface::GetBoundaryLimit()", "GeomSolids1002",
                 JustWarning, message);
   }
}

//=====================================================================
//* SetBoundary -------------------------------------------------------

void G4VTwistSurface::SetBoundary(const G4int         &axiscode,
                             const G4ThreeVector &direction,
                             const G4ThreeVector &x0,
                             const G4int         &boundarytype)
{
   G4int code = (~sAxisMask) & axiscode;
   if ((code == (sAxis0 & sAxisMin)) ||
       (code == (sAxis0 & sAxisMax)) ||
       (code == (sAxis1 & sAxisMin)) ||
       (code == (sAxis1 & sAxisMax))) {

      G4int i;
      G4bool done = false;
      for (i=0; i<4; i++) {
         if (fBoundaries[i].IsEmpty()) {
            fBoundaries[i].SetFields(axiscode, direction,
                                     x0, boundarytype);
            done = true;
            break;
         }
      }

      if (!done) {
         G4Exception("G4VTwistSurface::SetBoundary()", "GeomSolids0003",
                      FatalException, "Number of boundary exceeding 4!");
      }
   } else {
      std::ostringstream message;
      message << "Invalid axis-code." << G4endl
              << "        axiscode = "
              << std::hex << axiscode << std::dec;
      G4Exception("G4VTwistSurface::SetBoundary()", "GeomSolids0003",
                  FatalException, message);
   }
}

//=====================================================================
//* GetFace -----------------------------------------------------------

G4int G4VTwistSurface::GetFace( G4int i, G4int j, G4int k,
                                G4int n, G4int iside ) 
{
  // this is the face mapping function
  // (i,j) -> face number

  if ( iside == 0 ) {
    return i * ( k - 1 ) + j ;
  }

  else if ( iside == 1 ) {
    return (k-1)*(k-1) + i*(k-1) + j ;
  }

  else if ( iside == 2 ) {
    return 2*(k-1)*(k-1) + i*(k-1) + j ;
  }

  else if ( iside == 3 ) {
    return 2*(k-1)*(k-1) + (n-1)*(k-1) + i*(k-1) + j ;
  }
  
  else if ( iside == 4 ) {
    return 2*(k-1)*(k-1) + 2*(n-1)*(k-1) + i*(k-1) + j ;
  }
  
  else if ( iside == 5 ) {
    return 2*(k-1)*(k-1) + 3*(n-1)*(k-1) + i*(k-1) + j ;
  }

  else {

    std::ostringstream message;
    message << "Not correct side number: "
            << GetName() << G4endl
            << "iside is " << iside << " but should be "
            << "0,1,2,3,4 or 5" << ".";
    G4Exception("G4TwistSurface::G4GetFace()", "GeomSolids0002",
                FatalException, message);


  }

  return -1 ;  // wrong face
}

//=====================================================================
//* GetNode -----------------------------------------------------------

G4int G4VTwistSurface::GetNode( G4int i, G4int j, G4int k,
                                G4int n, G4int iside ) 
{
  // this is the node mapping function
  // (i,j) -> node number
  // Depends on the side iside and the used meshing of the surface

  if ( iside == 0 ) {
    // lower endcap is kxk squared. 
    // n = k 
    return i * k + j ;
  }

  if ( iside == 1 ) {
    // upper endcap is kxk squared. Shift by k*k
    // n = k 
    return  k*k + i*k + j ;
  }

  else if ( iside == 2 ) {
    // front side.
    if      ( i == 0 )     {   return       j ;  }
    else if ( i == n-1 )   {   return k*k + j ;  } 
    else                   {   return 2*k*k + 4*(i-1)*(k-1) + j ; }
  }

  else if ( iside == 3 ) {
    // right side
    if      ( i == 0 )     {   return       (j+1)*k - 1 ; } 
    else if ( i == n-1 )   {   return k*k + (j+1)*k - 1 ; }
    else                   {   return 2*k*k + 4*(i-1)*(k-1) + (k-1) + j ; }
  }
  else if ( iside == 4 ) {
    // back side.
    if      ( i == 0 )     {   return   k*k - 1 - j ; }               // reversed order
    else if ( i == n-1 )   {   return 2*k*k - 1 - j ; }               // reversed order 
    else                   {   return 2*k*k + 4*(i-1)*(k-1) + 2*(k-1) + j ; // normal order
    }
  }
  else if ( iside == 5 ) {
    // left side 
    if      ( i == 0 )     {   return k*k   - (j+1)*k ; }             // reversed order
    else if ( i == n-1)    {   return 2*k*k - (j+1)*k ; }             // reverded order
    else {
      if ( j == k-1 )      {   return 2*k*k + 4*(i-1)*(k-1) ; }       // special case
      else                 {   return 2*k*k + 4*(i-1)*(k-1) + 3*(k-1) + j ; }  // normal order
    }  
  }

  else {

    std::ostringstream message;
    message << "Not correct side number: "
            << GetName() << G4endl
            << "iside is " << iside << " but should be "
            << "0,1,2,3,4 or 5" << ".";
    G4Exception("G4TwistSurface::G4GetNode()", "GeomSolids0002",
                FatalException, message);
  } 
  return -1 ;  // wrong node 
} 

//=====================================================================
//* GetEdgeVisiblility ------------------------------------------------

G4int G4VTwistSurface::GetEdgeVisibility( G4int i, G4int j, G4int k, G4int n, G4int number, G4int orientation) 
{

  // clockwise filling         -> positive orientation
  // counter clockwise filling -> negative orientation 

  //
  //   d    C    c
  //     +------+ 
  //     |      |
  //     |      |
  //     |      |
  //   D |      |B
  //     |      |
  //     |      |
  //     |      |
  //     +------+
  //    a   A    b
  //    
  //  a = +--+    A = ---+
  //  b = --++    B = --+-
  //  c = -++-    C = -+--
  //  d = ++--    D = +---


  // check first invisible faces

  if ( ( i>0 && i<n-2 ) && ( j>0 && j<k-2 ) ) {
    return -1 ;   // always invisible, signs:   ----
  }
  
  // change first the vertex number (depends on the orientation)
  // 0,1,2,3  -> 3,2,1,0
  if ( orientation < 0 ) { number = ( 3 - number ) ;  }
  
  // check true edges
  if ( ( j>=1 && j<=k-3 ) ) {

    if ( i == 0 )  {        // signs (A):  ---+  
      return ( number == 3 ) ? 1 : -1 ;
    }
    
    else if ( i == n-2 ) {  // signs (C):  -+--
      return  ( number == 1 ) ? 1 : -1 ; 
    }
    
    else {
      std::ostringstream message;
      message << "Not correct face number: " << GetName() << " !";
      G4Exception("G4TwistSurface::G4GetEdgeVisibility()",
                  "GeomSolids0003", FatalException, message);
    }
  }
  
  if ( ( i>=1 && i<=n-3 ) ) {

    if ( j == 0 )  {        // signs (D):  +---
      return ( number == 0 ) ? 1 : -1 ;
    }

    else if ( j == k-2 ) {  // signs (B):  --+-
      return  ( number == 2 ) ? 1 : -1 ; 
    }

    else {
      std::ostringstream message;
      message << "Not correct face number: " << GetName() << " !";
      G4Exception("G4TwistSurface::G4GetEdgeVisibility()",
                  "GeomSolids0003", FatalException, message);
    }
  }
  
  // now the corners
  if ( i == 0 && j == 0 ) {          // signs (a) : +--+
    return ( number == 0 || number == 3 ) ? 1 : -1 ;
  }
  else if ( i == 0 && j == k-2 ) {   // signs (b) : --++
    return ( number == 2 || number == 3 ) ? 1 : -1 ;
  }
  else if ( i == n-2 && j == k-2 ) { // signs (c) : -++-
    return ( number == 1 || number == 2 ) ? 1 : -1 ;
  }
  else if ( i == n-2 && j == 0 ) {   // signs (d) : ++--
    return ( number == 0 || number == 1 ) ? 1 : -1 ;
  }
  else {
    std::ostringstream message;
    message << "Not correct face number: " << GetName() << " !";
    G4Exception("G4TwistSurface::G4GetEdgeVisibility()",
                "GeomSolids0003", FatalException, message);
  }

  std::ostringstream message;
  message << "Not correct face number: " << GetName() << " !";
  G4Exception("G4TwistSurface::G4GetEdgeVisibility()", "GeomSolids0003",
              FatalException, message);

  return 0 ;
}


//=====================================================================
//* DebugPrint --------------------------------------------------------

void G4VTwistSurface::DebugPrint() const
{
   G4ThreeVector A = fRot * GetCorner(sC0Min1Min) + fTrans;
   G4ThreeVector B = fRot * GetCorner(sC0Max1Min) + fTrans;
   G4ThreeVector C = fRot * GetCorner(sC0Max1Max) + fTrans;
   G4ThreeVector D = fRot * GetCorner(sC0Min1Max) + fTrans;
  
   G4cout << "/* G4VTwistSurface::DebugPrint():-------------------------------"
          << G4endl;
   G4cout << "/* Name = " << fName << G4endl;
   G4cout << "/* Axis = " << std::hex << fAxis[0] << " "
          << std::hex << fAxis[1] 
          << " (0,1,2,3,5 = kXAxis,kYAxis,kZAxis,kRho,kPhi)"
          << std::dec << G4endl;
   G4cout << "/* BoundaryLimit(in local) fAxis0(min, max) = ("<<fAxisMin[0] 
          << ", " << fAxisMax[0] << ")" << G4endl;
   G4cout << "/* BoundaryLimit(in local) fAxis1(min, max) = ("<<fAxisMin[1] 
          << ", " << fAxisMax[1] << ")" << G4endl;
   G4cout << "/* Cornar point sC0Min1Min = " << A << G4endl;
   G4cout << "/* Cornar point sC0Max1Min = " << B << G4endl;
   G4cout << "/* Cornar point sC0Max1Max = " << C << G4endl;
   G4cout << "/* Cornar point sC0Min1Max = " << D << G4endl;
   G4cout << "/*---------------------------------------------------------"
          << G4endl;
}

//=====================================================================
// G4VTwistSurface::CurrentStatus class
//=====================================================================

//=====================================================================
//* CurrentStatus::CurrentStatus --------------------------------------

G4VTwistSurface::CurrentStatus::CurrentStatus() 
{
  for (size_t i=0; i<G4VSURFACENXX; i++)
  {
    fDistance[i] = kInfinity;
    fAreacode[i] = sOutside;
    fIsValid[i]  = false;
    fXX[i].set(kInfinity, kInfinity, kInfinity);
  }
  fNXX   = 0;
  fLastp.set(kInfinity, kInfinity, kInfinity);
  fLastv.set(kInfinity, kInfinity, kInfinity);
  fLastValidate = kUninitialized;
  fDone = false;
}

//=====================================================================
//* CurrentStatus::~CurrentStatus -------------------------------------

G4VTwistSurface::CurrentStatus::~CurrentStatus() 
{
}

//=====================================================================
//* CurrentStatus::SetCurrentStatus -----------------------------------

void
G4VTwistSurface::CurrentStatus::SetCurrentStatus(G4int                i, 
                                            G4ThreeVector       &xx, 
                                            G4double            &dist, 
                                            G4int               &areacode, 
                                            G4bool              &isvalid,
                                            G4int                nxx,
                                            EValidate            validate,
                                      const G4ThreeVector *p, 
                                      const G4ThreeVector *v)
{
  fDistance[i]  = dist;
  fAreacode[i]  = areacode;
  fIsValid[i]   = isvalid;
  fXX[i]        = xx;
  fNXX          = nxx;
  fLastValidate = validate;
  if (p)
  {
    fLastp = *p;
  }
  else
  {
    G4Exception("G4VTwistSurface::CurrentStatus::SetCurrentStatus()",
                "GeomSolids0003", FatalException, "SetCurrentStatus: p = 0!");
  }
  if (v) 
  {
    fLastv = *v;
  }
  else
  {
    fLastv.set(kInfinity, kInfinity, kInfinity);
  }
  fDone = true;
}

//=====================================================================
//* CurrentStatus::ResetfDone -----------------------------------------

void
G4VTwistSurface::CurrentStatus::ResetfDone(EValidate     validate,
                                const G4ThreeVector *p, 
                                const G4ThreeVector *v)

{
  if (validate == fLastValidate && p && *p == fLastp)
  {
     if (!v || (*v == fLastv)) return;
  }         
  G4ThreeVector xx(kInfinity, kInfinity, kInfinity);
  for (size_t i=0; i<G4VSURFACENXX; i++)
  {
    fDistance[i] = kInfinity;
    fAreacode[i] = sOutside;
    fIsValid[i]  = false;
    fXX[i] = xx;   // bug in old code ( was fXX[i] =  xx[i]  )
  }
  fNXX   = 0;
  fLastp.set(kInfinity, kInfinity, kInfinity);
  fLastv.set(kInfinity, kInfinity, kInfinity);
  fLastValidate = kUninitialized;
  fDone = false;
}

//=====================================================================
//* CurrentStatus::DebugPrint -----------------------------------------

void
G4VTwistSurface::CurrentStatus::DebugPrint() const
{
  G4cout << "CurrentStatus::Dist0,1= " << fDistance[0] 
         << " " << fDistance[1] << " areacode = " << fAreacode[0] 
         << " " << fAreacode[1] << G4endl;
}

//=====================================================================
// G4VTwistSurface::Boundary class
//=====================================================================

//=====================================================================
//* Boundary::Boundary ------------------------------------------------

G4VTwistSurface::Boundary::Boundary()
 : fBoundaryAcode(-1), fBoundaryType(0)
{
}

//=====================================================================
//* Boundary::~Boundary -----------------------------------------------

G4VTwistSurface::Boundary::~Boundary()
{
}

//=====================================================================
//* Boundary::SetFields -----------------------------------------------

void
G4VTwistSurface::Boundary::SetFields(const G4int         &areacode, 
                                const G4ThreeVector &d, 
                                const G4ThreeVector &x0, 
                                const G4int         &boundarytype)
{
  fBoundaryAcode     = areacode;
  fBoundaryDirection = d;
  fBoundaryX0        = x0;
  fBoundaryType      = boundarytype;
}

//=====================================================================
//* Boundary::IsEmpty -------------------------------------------------

G4bool G4VTwistSurface::Boundary::IsEmpty() const 
{
  if (fBoundaryAcode == -1) return true;
  return false;
}

//=====================================================================
//* Boundary::GetBoundaryParameters -----------------------------------

G4bool
G4VTwistSurface::Boundary::GetBoundaryParameters(const G4int         &areacode, 
                                                  G4ThreeVector &d,
                                                  G4ThreeVector &x0, 
                                                  G4int  &boundarytype) const 
{  
  // areacode must be one of them:
  // sAxis0 & sAxisMin, sAxis0 & sAxisMax,
  // sAxis1 & sAxisMin, sAxis1 & sAxisMax
  if ((areacode & sAxis0) && (areacode & sAxis1))
  {
    std::ostringstream message;
    message << "Located in the corner area." << G4endl
            << "        This function returns a direction vector of "
            << "a boundary line." << G4endl
            << "        areacode = " << areacode;
    G4Exception("G4VTwistSurface::Boundary::GetBoundaryParameters()",
                "GeomSolids0003", FatalException, message);
  } 
  if ((areacode & sSizeMask) != (fBoundaryAcode & sSizeMask))
  {
    return false;
  }
  d  = fBoundaryDirection;
  x0 = fBoundaryX0;
  boundarytype = fBoundaryType;
  return true;
}
