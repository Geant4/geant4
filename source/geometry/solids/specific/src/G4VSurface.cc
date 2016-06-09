//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VSurface.cc,v 1.11 2004/12/17 16:34:59 link Exp $
// GEANT4 tag $Name: geant4-07-00-patch-01 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4VSurface.cc
//
// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------

#include "G4VSurface.hh"
#include <iomanip>

const G4int  G4VSurface::sOutside        = 0x00000000;
const G4int  G4VSurface::sInside         = 0x10000000;
const G4int  G4VSurface::sBoundary       = 0x20000000;
const G4int  G4VSurface::sCorner         = 0x40000000;
const G4int  G4VSurface::sC0Min1Min      = 0x40000101; 
const G4int  G4VSurface::sC0Max1Min      = 0x40000201;
const G4int  G4VSurface::sC0Max1Max      = 0x40000202; 
const G4int  G4VSurface::sC0Min1Max      = 0x40000102; 
const G4int  G4VSurface::sAxisMin        = 0x00000101; 
const G4int  G4VSurface::sAxisMax        = 0x00000202; 
const G4int  G4VSurface::sAxisX          = 0x00000404;
const G4int  G4VSurface::sAxisY          = 0x00000808;
const G4int  G4VSurface::sAxisZ          = 0x00000C0C;
const G4int  G4VSurface::sAxisRho        = 0x00001010;
const G4int  G4VSurface::sAxisPhi        = 0x00001414;

// mask
const G4int  G4VSurface::sAxis0          = 0x0000FF00;
const G4int  G4VSurface::sAxis1          = 0x000000FF;
const G4int  G4VSurface::sSizeMask       = 0x00000303;
const G4int  G4VSurface::sAxisMask       = 0x0000FCFC;
const G4int  G4VSurface::sAreaMask       = 0XF0000000;

//=====================================================================
//* constructors ------------------------------------------------------

G4VSurface::G4VSurface(const G4String &name)
  : fName(name)
{

   fAxis[0]    = kUndefined;
   fAxis[1]    = kUndefined;
   fAxisMin[0] = kInfinity;
   fAxisMin[1] = kInfinity;
   fAxisMax[0] = kInfinity;
   fAxisMax[1] = kInfinity;
   fHandedness = 1;

   G4int i;
   for (i=0; i<4; i++) {
      fCorners[i].set(kInfinity, kInfinity, kInfinity);
      fNeighbours[i] = 0;
   }

   fCurrentNormal.p.set(kInfinity, kInfinity, kInfinity);
   
   fAmIOnLeftSide.me.set(kInfinity, kInfinity, kInfinity);
   fAmIOnLeftSide.vec.set(kInfinity, kInfinity, kInfinity);
}

G4VSurface::G4VSurface(const G4String         &name,
                       const G4RotationMatrix &rot,
                       const G4ThreeVector    &tlate,
                             G4int             handedness,
                       const EAxis             axis0 ,
                       const EAxis             axis1 ,
                             G4double          axis0min,
                             G4double          axis1min,
                             G4double          axis0max,
                             G4double          axis1max )
   : fName(name)
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

   G4int i;
   for (i=0; i<4; i++) {
      fCorners[i].set(kInfinity, kInfinity, kInfinity);
      fNeighbours[i] = 0;
   }

   fCurrentNormal.p.set(kInfinity, kInfinity, kInfinity);
   
   fAmIOnLeftSide.me.set(kInfinity, kInfinity, kInfinity);
   fAmIOnLeftSide.vec.set(kInfinity, kInfinity, kInfinity);
}

//=====================================================================
//* destructor --------------------------------------------------------

G4VSurface::~G4VSurface()
{
}

//=====================================================================
//* AmIOnLeftSide -----------------------------------------------------

G4int G4VSurface::AmIOnLeftSide(const G4ThreeVector &me, 
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

   static G4RotationMatrix unitrot;  // unit matrix
   static const G4RotationMatrix rottol    = unitrot.rotateZ(0.5*kAngTolerance);
   static const G4RotationMatrix invrottol = unitrot.rotateZ(-1.*kAngTolerance);

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

#ifdef G4SPECSDEBUG
   G4cout << "         === G4VSurface::AmIOnLeftSide() =============="
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

G4double G4VSurface::DistanceToBoundary(G4int areacode,
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
      G4cerr << "ERROR - G4VSurface::DistanceToBoundary()" << G4endl
             << "        Point is in the corner area. This function returns"
             << G4endl
             << "        a direction vector of a boundary line." << G4endl
             << "        areacode = " << areacode << G4endl;
      G4Exception("G4VSurface::DistanceToBoundary()", "InvalidSetup",
                  FatalException, "Point is in the corner area.");
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
      G4cerr << "ERROR - G4VSurface::DistanceToBoundary()" << G4endl
             << "        areacode = " << areacode << G4endl;
      G4Exception("G4VSurface::DistanceToBoundary()", "InvalidSetup",
                  FatalException, "Bad areacode of boundary.");
   }
   return dist;
}

//=====================================================================
//* DistanceToIn ------------------------------------------------------

G4double G4VSurface::DistanceToIn(const G4ThreeVector &gp,
                                  const G4ThreeVector &gv,
                                        G4ThreeVector &gxxbest)
{
#ifdef G4SPECSDEBUG
   G4cout << " ~~~~~ G4VSurface::DistanceToIn(p,v) - Start ~~~~~" << G4endl;
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
   G4int         besti          = -1;  
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

#ifdef G4SPECSDEBUG
         G4cout << "   G4VSurface::DistanceToIn(p,v): "
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
            besti   = i;

#ifdef G4SPECSDEBUG
            G4cout << "   G4VSurface::DistanceToIn(p,v): "
                   << " areacode sInside name, distance = "
                   << fName <<  " "<< bestdistance << G4endl;
#endif 
         }

      //
      // else, the intersection is on boundary or corner.
      //

      } else {

         G4VSurface *neighbours[2];
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

#ifdef G4SPECSDEBUG
                  G4cout << "   G4VSurface:DistanceToIn(p,v): "
                         << " intersection "<< tmpgxx[k] << G4endl
                         << "   is inside of neighbour surface of " << fName 
                         << " . returning kInfinity." << G4endl;
                  G4cout << "~~~~~ G4VSurface::DistanceToIn(p,v) - return ~~~~"
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
                besti   = i;
#ifdef G4SPECSDEBUG
               G4cout << "   G4VSurface::DistanceToIn(p,v): "
                      << " areacode sBoundary & sBoundary distance = "
                      << fName  << " " << distance[i] << G4endl;
#endif 
            }
         }

      } // else end
   } // intersection loop end

   gxxbest = bestgxx;

#ifdef G4SPECSDEBUG
   if (besti < 0) {
      G4cout << "~~~~~ G4VSurface::DistanceToIn(p,v) - return ~~~~" << G4endl;
      G4cout << "      No intersections " << G4endl; 
      G4cout << "      Name : " << fName << G4endl; 
      G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
   } else {
      G4cout << "~~~~~ G4VSurface::DistanceToIn(p,v) : return ~~~~" << G4endl;
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

G4double G4VSurface::DistanceToOut(const G4ThreeVector &gp,
                                   const G4ThreeVector &gv,
                                         G4ThreeVector &gxxbest)
{
#ifdef G4SPECSDEBUG
   G4cout << "~~~~~ G4VSurface::DistanceToOut(p,v) - Start ~~~~" << G4endl;
   G4cout << "      Name : " << fName << G4endl;
   G4cout << "      gp   : " << gp << G4endl;
   G4cout << "      gv   : " <<  gv << G4endl;
   G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
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

   G4int         nxx;
   G4double      bestdistance   = kInfinity;
   G4int         besti          = -1;

   nxx = DistanceToSurface(gp, gv, gxx, distance, areacode,
                           isvalid, kValidateWithTol);
   G4int i;
   for (i=0; i<nxx; i++) {
      if (!(isvalid[i])) {
         continue;
      }

      G4ThreeVector normal = GetNormal(gxx[i], true);
      if (normal * gv <= 0) {
         // particle goes toword inside of solid, return kInfinity
#ifdef G4SPECSDEBUG
          G4cout << "   G4VSurface::DistanceToOut(p,v): normal*gv < 0, normal " 
                 << fName << " " << normal 
                 << G4endl;
#endif 
      } else {
         // gxx[i] is accepted.
         if (distance[i] < bestdistance) {
            bestdistance = distance[i];
            gxxbest = gxx[i];
            besti   = i;
         }
      } 
   }

#ifdef G4SPECSDEBUG
   if (besti < 0) {
      G4cout << "~~~~~ G4VSurface::DistanceToOut(p,v) - return ~~~" << G4endl;
      G4cout << "      No intersections   " << G4endl; 
      G4cout << "      Name     : " << fName << G4endl; 
      G4cout << "      bestdist : " << bestdistance << G4endl;
      G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
   } else {
      G4cout << "~~~~~ G4VSurface::DistanceToOut(p,v) : return ~~~" << G4endl;
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

G4double G4VSurface::DistanceTo(const G4ThreeVector &gp,
                                      G4ThreeVector &gxxbest)
{
#ifdef G4SPECSDEBUG
   G4cout << "~~~~~ G4VSurface::DistanceTo(p) - Start ~~~~~~~~~" << G4endl;
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

   G4int nxx;

   nxx = DistanceToSurface(gp, gxx, distance, areacode);
   gxxbest = gxx[0];

#ifdef G4SPECSDEBUG
   G4cout << "~~~~~ G4VSurface::DistanceTo(p) - return ~~~~~~~~" << G4endl;
   G4cout << "      Name     : " << fName << G4endl; 
   G4cout << "      gxx      : " << gxxbest << G4endl; 
   G4cout << "      bestdist : " << distance[0] << G4endl;
   G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
#endif

   return distance[0];
}

//=====================================================================
//* IsSameBoundary ----------------------------------------------------

G4bool G4VSurface::IsSameBoundary(G4VSurface *surface1, G4int areacode1,
                                  G4VSurface *surface2, G4int areacode2 ) const
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

void G4VSurface::GetBoundaryParameters(const G4int         &areacode,
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

   G4cerr << "ERROR - G4VSurface::GetBoundaryParameters()" << G4endl
          << "        Boundary at areacode " << std::hex << areacode
          << std::dec << G4endl
          << "        is not be registered." << G4endl;
   G4Exception("G4VSurface::GetBoundaryParameters()", "InvalidSetup",
               FatalException, "Not registered boundary.");
}

//=====================================================================
//* GetBoundaryAtPZ ---------------------------------------------------

G4ThreeVector G4VSurface::GetBoundaryAtPZ(G4int areacode,
                                          const G4ThreeVector &p) const
{
   // areacode must be one of them:
   // sAxis0 & sAxisMin, sAxis0 & sAxisMax,
   // sAxis1 & sAxisMin, sAxis1 & sAxisMax.

   if (areacode & sAxis0 && areacode & sAxis1) {
     G4cerr << "ERROR - G4VSurface::GetBoundaryAtPZ()" << G4endl
            << "        Point is in the corner area. This function returns"
            << G4endl
            << "        a direction vector of a boundary line." << G4endl
            << "        areacode = " << areacode << G4endl;
     G4Exception("G4VSurface::GetBoundaryAtPZ()", "InvalidCondition",
                 FatalException, "Point is in the corner area.");
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
     G4cerr << "ERROR - G4VSurface::GetBoundaryAtPZ()" << G4endl
            << "        Boundary at areacode " << areacode << G4endl
            << "        is not be registered." << G4endl;
     G4Exception("G4VSurface::GetBoundaryAtPZ()", "InvalidSetup",
                 FatalException, "Not registered boundary.");
   }

   if (((boundarytype & sAxisPhi) == sAxisPhi) ||
       ((boundarytype & sAxisRho) == sAxisRho)) {
     G4cerr << "ERROR - G4VSurface::GetBoundaryAtPZ()" << G4endl
            << "        Boundary at areacode " << areacode << G4endl
            << "        is not a z-depended line." << G4endl;
     G4Exception("G4VSurface::GetBoundaryAtPZ()", "InvalidSetup",
                 FatalException, "Not a z-depended line boundary.");
   }
   return ((p.z() - x0.z()) / d.z()) * d + x0;
}

//=====================================================================
//* SetCorner ---------------------------------------------------------

void G4VSurface::SetCorner(G4int areacode, G4double x, G4double y, G4double z)
{
   if ((areacode & sCorner) != sCorner){
     G4cerr << "ERROR - G4VSurface::SetCorner()" << G4endl
            << "        areacode " << areacode << G4endl;
     G4Exception("G4VSurface::SetCorner()", "InvalidSetup",
                 FatalException, "Area code must represents corner.");
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

void G4VSurface::GetBoundaryAxis(G4int areacode, EAxis axis[]) const
{
   if ((areacode & sBoundary) != sBoundary) {
     G4Exception("G4VSurface::GetBoundaryAxis()", "InvalidCondition",
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
           G4cerr << "ERROR - G4VSurface::GetBoundaryAxis()" << G4endl
                  << "        areacode " << areacode << G4endl;
           G4Exception("G4VSurface::GetBoundaryAxis()", "InvalidSetup",
                       FatalException, "Not supported areacode.");
         }
      }
   }
}

//=====================================================================
//* SetBoundaryLimit --------------------------------------------------

void G4VSurface::GetBoundaryLimit(G4int areacode, G4double limit[]) const
{
   if (areacode & sCorner) {
      if (areacode & sC0Min1Max) {
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
     G4cerr << "WARNING - G4VSurface::GetBoundaryAxis()" << G4endl
            << "          areacode " << areacode << G4endl;
     G4Exception("G4VSurface::GetBoundaryLimit()", "InvalidCondition",
                 JustWarning, "Not located on a boundary!");
   }
}

//=====================================================================
//* SetBoundary -------------------------------------------------------

void G4VSurface::SetBoundary(const G4int         &axiscode,
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
         G4Exception("G4VSurface::SetBoundary()", "InvalidCondition",
                      FatalException, "Number of boundary exceeding 4!");
      }
   } else {
      G4cerr << "ERROR - G4VSurface::SetBoundary()" << G4endl
             << "        invalid axiscode. axiscode = "
             << std::hex << axiscode << std::dec << G4endl;
      G4Exception("G4VSurface::SetBoundary()", "InvalidCondition",
                  FatalException, "Invalid axis-code.");
   }
}

//=====================================================================
//* DebugPrint --------------------------------------------------------

void G4VSurface::DebugPrint() const
{
   G4ThreeVector A = fRot * GetCorner(sC0Min1Min) + fTrans;
   G4ThreeVector B = fRot * GetCorner(sC0Max1Min) + fTrans;
   G4ThreeVector C = fRot * GetCorner(sC0Max1Max) + fTrans;
   G4ThreeVector D = fRot * GetCorner(sC0Min1Max) + fTrans;
  
   G4cout << "/* G4VSurface::DebugPrint():-------------------------------"
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
// G4VSurface::CurrentStatus class
//=====================================================================

//=====================================================================
//* CurrentStatus::CurrentStatus --------------------------------------

G4VSurface::CurrentStatus::CurrentStatus() 
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

G4VSurface::CurrentStatus::~CurrentStatus() 
{
}

//=====================================================================
//* CurrentStatus::SetCurrentStatus -----------------------------------

void
G4VSurface::CurrentStatus::SetCurrentStatus(G4int                i, 
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
    G4Exception("G4VSurface::CurrentStatus::CurrentStatus()",
                "InvalidCondition", FatalException,
                "SetCurrentStatus: p = 0!");
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
G4VSurface::CurrentStatus::ResetfDone(EValidate     validate,
                                const G4ThreeVector *p, 
                                const G4ThreeVector *v)

{
  if (validate == fLastValidate && p && *p == fLastp)
  {
     if (!v || (v && *v == fLastv)) return;
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
G4VSurface::CurrentStatus::DebugPrint() const
{
  G4cout << "CurrentStatus::Dist0,1= " << fDistance[0] 
         << " " << fDistance[1] << " areacode = " << fAreacode[0] 
         << " " << fAreacode[1] << G4endl;
}

//=====================================================================
// G4VSurface::Boundary class
//=====================================================================

//=====================================================================
//* Boundary::Boundary ------------------------------------------------

G4VSurface::Boundary::Boundary()
 : fBoundaryAcode(-1), fBoundaryType(0)
{
}

//=====================================================================
//* Boundary::~Boundary -----------------------------------------------

G4VSurface::Boundary::~Boundary()
{
}

//=====================================================================
//* Boundary::SetFields -----------------------------------------------

void
G4VSurface::Boundary::SetFields(const G4int         &areacode, 
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

G4bool G4VSurface::Boundary::IsEmpty() const 
{
  if (fBoundaryAcode == -1) return true;
  return false;
}

//=====================================================================
//* Boundary::GetBoundaryParameters -----------------------------------

G4bool
G4VSurface::Boundary::GetBoundaryParameters(const G4int         &areacode, 
                                                  G4ThreeVector &d,
                                                  G4ThreeVector &x0, 
                                                  G4int  &boundarytype) const 
{  
  // areacode must be one of them:
  // sAxis0 & sAxisMin, sAxis0 & sAxisMax,
  // sAxis1 & sAxisMin, sAxis1 & sAxisMax
  if ((areacode & sAxis0) && (areacode & sAxis1))
  {
    G4cerr << "ERROR - G4VSurface::Boundary::GetBoundaryParameters()"
           << G4endl
           << "        Located in the corner area. This function"
           << "        returns a direction vector of a boundary line."
           << "        areacode = " << areacode << G4endl;
    G4Exception("G4VSurface::Boundary::GetBoundaryParameters()",
                "InvalidCondition", FatalException,
                "Located in the corner area.");
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
