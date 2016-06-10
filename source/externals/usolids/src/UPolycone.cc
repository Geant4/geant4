//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UPolycone
//
// 19.04.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include "UPolycone.hh"

#include "UEnclosingCylinder.hh"
#include "UReduciblePolygon.hh"

#include "UTubs.hh"
#include "UCons.hh"
#include "UTransform3D.hh"

using namespace std;

UPolycone::UPolycone(const std::string& name,
                     double phiStart,
                     double phiTotal,
                     int numZPlanes,
                     const double zPlane[],
                     const double rInner[],
                     const double rOuter[])
  : VUSolid(name)  //, fNumSides(0)
{
  fCubicVolume = 0;
  fSurfaceArea = 0;
  Init(phiStart, phiTotal, numZPlanes, zPlane, rInner, rOuter);

}

UPolycone::UPolycone(const std::string& name,
                     double phiStart,
                     double phiTotal,
                     int    numRZ,
                     const double r[],
                     const double z[])
  : VUSolid(name)
{
  UReduciblePolygon* rz = new UReduciblePolygon(r, z, numRZ);

  // Create( phiStart, phiTotal, rz );

  // Set original_parameters struct for consistency
  //

  bool convertible = SetOriginalParameters(rz);

  if (!convertible)
  {
    std::ostringstream message;
    message << "Polycone " << GetName() << "cannot be converted" << std::endl
            << "to Polycone with (Rmin,Rmaz,Z) parameters! Use GenericPolycone" ;
    UUtils::Exception("UPolycone::UPolycone()", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }
  else
  {
    std::cout << "INFO: Converting polycone " << GetName() << std::endl
              << "to optimized polycone with (Rmin,Rmaz,Z) parameters !"
              << std::endl;
    double* Z, *R1, *R2;
    int num = fOriginalParameters->fNumZPlanes;
    Z = new double[num];
    R1 = new double[num];
    R2 = new double[num];
    for (int i = 0; i < num; i++)
    {
      Z[i] = fOriginalParameters->fZValues[i];
      R1[i] = fOriginalParameters->Rmin[i];
      R2[i] = fOriginalParameters->Rmax[i];
    }

    delete(fOriginalParameters);
    Init(phiStart, phiTotal, num, Z, R1, R2);
    delete [] R1;
    delete [] Z;
    delete [] R2;
  }

  delete rz;
}


//
// Constructor (GEANT3 style parameters)
//
void UPolycone::Init(double phiStart,
                     double phiTotal,
                     int numZPlanes,
                     const double zPlane[],
                     const double rInner[],
                     const double rOuter[])
{
  //Convertion for angles
 
  if (phiTotal <= 0 || phiTotal > UUtils::kTwoPi-1E-10)
   {
     phiIsOpen=false;
     startPhi = 0;
     endPhi = UUtils::kTwoPi;
   }
   else
   {
     //
     // Convert phi into our convention
     //
     phiIsOpen=true;
     startPhi = phiStart;
     while( startPhi < 0 ) startPhi += UUtils::kTwoPi;
     
     endPhi = phiStart+phiTotal;
     while( endPhi < startPhi ) endPhi += UUtils::kTwoPi;

     
  }
 
  // Set Parameters
  fOriginalParameters = new UPolyconeHistorical();
  fOriginalParameters->fStartAngle = startPhi;
  fOriginalParameters->fOpeningAngle = endPhi-startPhi;
  fOriginalParameters->fNumZPlanes = numZPlanes;
  fOriginalParameters->fZValues.resize(numZPlanes);
  fOriginalParameters->Rmin.resize(numZPlanes);
  fOriginalParameters->Rmax.resize(numZPlanes);

  // Calculate RMax of Polycone in order to determine convexity of sections
  //
  double RMaxextent=rOuter[0];
  for (int j=1; j < numZPlanes; j++)
  {
    if (rOuter[j] > RMaxextent) RMaxextent=rOuter[j];
    if (rInner[j]>rOuter[j])
    {
      std::ostringstream message;
      message << "Cannot create Polycone with rInner > rOuter for the same Z"
              << std::endl
              << "        rInner > rOuter for the same Z !" << std::endl
              << "        rMin[" << j << "] = " << rInner[j]
              << " -- rMax[" << j << "] = " << rOuter[j];
       UUtils::Exception("UPolycone::UPolycone()", "GeomSolids0002",
                         UFatalErrorInArguments, 1, message.str().c_str());
    }
  }
  //
  double prevZ = zPlane[0], prevRmax = 0, prevRmin = 0;
  int dirZ = 1;
  if (zPlane[1] < zPlane[0])dirZ = -1;
//  int curSolid = 0;

  int i;
  for (i = 0; i < numZPlanes; i++)
  {
    if ((i < numZPlanes - 1) && (zPlane[i] == zPlane[i + 1]))
    {
      if ((rInner[i]  > rOuter[i + 1])
          || (rInner[i + 1] > rOuter[i]))
      {

        std::ostringstream message;
        message << "Cannot create a Polycone with no contiguous segments."
                << std::endl
                << "Segments are not contiguous !" << std::endl
                << " rMin[" << i << "] = " << rInner[i]
                << " -- rMax[" << i + 1 << "] = " << rOuter[i + 1] << std::endl
                << " rMin[" << i + 1 << "] = " << rInner[i + 1]
                << " -- rMax[" << i << "] = " << rOuter[i];
        UUtils::Exception("UPolycone::UPolycone()", "GeomSolids0002",
                          UFatalErrorInArguments, 1, message.str().c_str());
      }
    }



    double rMin = rInner[i];
    double rMax = rOuter[i];
    double z = zPlane[i];

    if (i > 0)
    {
      if (((z > prevZ)&&(dirZ>0))||((z < prevZ)&&(dirZ<0)))
      {
        if (dirZ*(z-prevZ)< 0)
        {
          std::ostringstream message;
          message << "Cannot create a Polycone with different Z directions.Use GenericPolycone."
                  << std::endl
                  << "ZPlane is changing direction  !" << std::endl
                  << "  zPlane[0] = " << zPlane[0]
                  << " -- zPlane[1] = " << zPlane[1] << std::endl
                  << "  zPlane[" << i - 1 << "] = " << zPlane[i - 1]
                  << " -- rPlane[" << i << "] = " << zPlane[i];
          UUtils::Exception("UPolycone::UPolycone()", "GeomSolids0002",
                            UFatalErrorInArguments, 1, message.str().c_str());



        }
        VUSolid* solid;
        double dz =  (z - prevZ)*dirZ / 2;

        bool tubular = (rMin == prevRmin && prevRmax == rMax);

//        if (fNumSides == 0)
        {
          if (tubular)
          {
            solid = new UTubs("", rMin, rMax, dz, phiStart, phiTotal);
          }
          else
          {
            solid = new UCons("", prevRmin, prevRmax, rMin, rMax, dz, phiStart, phiTotal);
          }
        }
//        else
//        {
//          solid = new UHedra("", prevRmin, prevRmax, rMin, rMax, dz, phiStart, phiTotal, fNumSides);
//        }

        fZs.push_back(z);

        int zi = fZs.size() - 1;
        double shift = fZs[zi - 1] + 0.5 * (fZs[zi] - fZs[zi - 1]);

        UPolyconeSection section;
        section.shift = shift;
        section.tubular = tubular;
        section.solid = solid;
        if(tubular)
        {
          if (rMax < RMaxextent) { section.convex = false;}
          else { section.convex = true;}
        }
        else
        {
          if ((rMax<prevRmax)||(rMax < RMaxextent)||(prevRmax < RMaxextent))
            { section.convex = false;}
          else
            { section.convex = true;}
        }
        fSections.push_back(section);
      }
      else
      {
        ;// i = i;
      }
    }
    else fZs.push_back(z);

    fOriginalParameters->fZValues[i] = zPlane[i];
    fOriginalParameters->Rmin[i] = rInner[i];
    fOriginalParameters->Rmax[i] = rOuter[i];

    prevZ = z;
    prevRmin = rMin;
    prevRmax = rMax;
  }

  fMaxSection = fZs.size() - 2;

  //
  // Build RZ polygon using special PCON/PGON GEANT3 constructor
  //
  UReduciblePolygon* rz = new UReduciblePolygon(rInner, rOuter, zPlane, numZPlanes);

  double mxy = rz->Amax();
//  double alfa = UUtils::kPi / fNumSides;

  double r = rz->Amax();
//
// Perform checks of rz values
//
  if (rz->Amin() < 0.0)
  {
     std::ostringstream message;
     message << "Illegal input parameters - " << GetName() << std::endl
             << "        All R values must be >= 0 !";
     UUtils::Exception("UPolycone::Init()", "GeomSolids0002",
                       UFatalErrorInArguments,1, message.str().c_str());
  }
     

  mxy += fgTolerance;

  fBox.Set(mxy, mxy, (rz->Bmax() - rz->Bmin()) / 2);

  //
  // Make enclosingCylinder
  //

  enclosingCylinder = new UEnclosingCylinder(r, rz->Bmax(), rz->Bmin(), phiIsOpen, phiStart, phiTotal);

  delete rz;
}


/*
//
// Constructor (generic parameters)
//
UPolycone3::UPolycone3( const std::string& name,
                              double phiStart,
                              double phiTotal,
                              int   numRZ,
                        const double r[],
                        const double z[]   )
  : VUSolid( name )
{
  UReduciblePolygon *rz = new UReduciblePolygon( r, z, numRZ );

  box.Set(rz->Amax(), rz->Amax(), (rz->Bmax() - rz->Bmin()) /2);

  // Set fOriginalParameters struct for consistency
  //
  SetOriginalParameters();

  delete rz;
}
*/


//
// Destructor
//
UPolycone::~UPolycone()
{
  //delete [] corners;
  //delete fOriginalParameters;
}



//
// Stream object contents to an output stream
//
std::ostream& UPolycone::StreamInfo(std::ostream& os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "                *** Dump for solid - " << GetName() << " ***\n"
     << "                ===================================================\n"
     << " Solid type: UPolycone\n"
     << " Parameters: \n"
     << "  starting phi angle : " << startPhi / (UUtils::kPi / 180.0) << " degrees \n"
     << "  ending phi angle   : " << endPhi / (UUtils::kPi / 180.0) << " degrees \n";
  int i = 0;
  int numPlanes = fOriginalParameters->fNumZPlanes;
    os << "  number of Z planes: " << numPlanes << "\n"
       << "            Z values: \n";
    for (i = 0; i < numPlanes; i++)
    {
      os << "    Z plane " << i << ": "
         << fOriginalParameters->fZValues[i] << "\n";
    }
    os << "  Tangent distances to inner surface (Rmin): \n";
    for (i = 0; i < numPlanes; i++)
    {
      os << "    Z plane " << i << ": "
         << fOriginalParameters->Rmin[i] << "\n";
    }
    os << "  Tangent distances to outer surface (Rmax): \n";
    for (i = 0; i < numPlanes; i++)
    {
      os << "    Z plane " << i << ": "
         << fOriginalParameters->Rmax[i] << "\n";
    }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

VUSolid::EnumInside UPolycone::InsideSection(int index, const UVector3& p) const
{
  const UPolyconeSection& section = fSections[index];
  UVector3 ps(p.x(), p.y(), p.z() - section.shift);

//  if (fNumSides) return section.solid->Inside(ps);

  double rMinPlus, rMaxPlus , rMinMinus, rMaxMinus;
  double dz;
  static const double halfTolerance = fgTolerance * 0.5;

  double r2 = p.x() * p.x() + p.y() * p.y();

  if (section.tubular)
  {
    UTubs* tubs = (UTubs*) section.solid;
     
    rMaxPlus = tubs->GetRMax()+ halfTolerance;
    rMinPlus = tubs->GetRMin() + halfTolerance;
    rMinMinus = tubs->GetRMin() - halfTolerance;
    rMaxMinus = tubs->GetRMax() - halfTolerance;
    dz = tubs->GetDz();//GetZHalfLength();
  }
  else
  {
    UCons* cons = (UCons*) section.solid;

    double rMax1 = cons->GetRmax1();
    double rMax2 = cons->GetRmax2();
    double rMin1 = cons->GetRmin1();
    double rMin2 = cons->GetRmin2();
    dz = cons->GetDz();
    double ratio = (ps.z() + dz) / (2 * dz);
    rMinPlus = rMin1 + (rMin2 - rMin1) * ratio + halfTolerance;
    rMaxPlus = rMax1 + (rMax2 - rMax1) * ratio + halfTolerance;
    rMinMinus = rMinPlus - 2*halfTolerance;
    rMaxMinus = rMaxPlus - 2*halfTolerance;
  }

  if (r2 < 1e-10)
  {
    EnumInside pos = eOutside;
    if ((!phiIsOpen) && (rMinMinus <= 0))
    {
      if(pos !=eSurface)  { pos = eInside; }
    }
    if (rMinMinus <= 0)  { pos = eSurface; }
    if ( ((pos == eInside) || (pos == eSurface))
      && ((ps.z() < -dz+halfTolerance) || (ps.z() > dz-halfTolerance)))
    {
      pos = eSurface;
    }
    return pos;    
  }
  
  if ( (r2 < rMinMinus*rMinMinus) || (r2 > rMaxPlus*rMaxPlus) )
  {
    return eOutside;
  }
  if ( (r2 < rMinPlus*rMinPlus) || (r2 > rMaxMinus*rMaxMinus) ) 
  {
    if(! phiIsOpen)  { return eSurface; }
  }

  if (! phiIsOpen)
  {
    if ( (ps.z() < -dz+halfTolerance) || (ps.z() > dz-halfTolerance) )
    {
      return eSurface;
    }
    return eInside;
  }
  double phi = std::atan2(p.y(), p.x()); // * UUtils::kTwoPi;
   
  if ( (phi < 0) || (endPhi > UUtils::kTwoPi) )  { phi += UUtils::kTwoPi; }

  double ddp = phi - startPhi;
  if (ddp < 0) ddp += UUtils::kTwoPi;
  
  if ( (phi <= endPhi+ frTolerance) && (phi>= startPhi-frTolerance) )
  { 
    if ( (ps.z() < -dz+halfTolerance) || (ps.z() > dz-halfTolerance) )
    {
      return eSurface;
    }
    if ( (r2 < rMinPlus*rMinPlus) || (r2 > rMaxMinus*rMaxMinus) )
    {
      return eSurface;
    }
    if (std::fabs(endPhi - phi) < frTolerance)
    {
      return eSurface;
    }
    if (std::fabs(startPhi - phi) < frTolerance)
    {
      return eSurface;
    }
    return eInside;
  }
  return eOutside;
}


VUSolid::EnumInside UPolycone::Inside(const UVector3& p) const
{
  double shift = fZs[0] + fBox.GetZHalfLength();
  UVector3 pb(p.x(), p.y(), p.z() - shift);

  if (fBox.Inside(pb) == eOutside)
  {
    return eOutside;
  }

  static const double htolerance = 0.5 * fgTolerance;
  int index = GetSection(p.z());

  EnumInside pos = InsideSection(index, p);
  if (pos == eInside)  { return eInside; }

  int nextSection;
  EnumInside nextPos;

  if (index > 0 && p.z()  - fZs[index] < htolerance)
  {
    nextSection = index - 1;
    nextPos = InsideSection(nextSection, p);
  }
  else if (index < fMaxSection && fZs[index + 1] - p.z() < htolerance)
  {
    nextSection = index + 1;
    nextPos = InsideSection(nextSection, p);
  }
  else
  {
    return pos;
  }

  if (nextPos == eInside)  { return eInside; }

  if (pos == eSurface && nextPos == eSurface)
  {
    UVector3 n, n2;
    NormalSection(index, p, n);
    NormalSection(nextSection, p, n2);
    if ((n +  n2).Mag2() < 1000 * frTolerance)  { return eInside; }
  }

  return (nextPos == eSurface || pos == eSurface) ? eSurface : eOutside;
}


double UPolycone::DistanceToIn(const UVector3& p,
                               const UVector3& v, double) const
{
  double shift = fZs[0] + fBox.GetZHalfLength();
  UVector3 pb(p.x(), p.y(), p.z() - shift);

  double idistance;

  idistance = fBox.DistanceToIn(pb, v);
    // using only box, this appears
    // to be faster than: idistance = enclosingCylinder->DistanceTo(pb, v);
  if (idistance >= UUtils::kInfinity)  { return idistance; }

  // this line can be here or not. not a big difference in performance
  // TODO: fix enclosingCylinder for polyhedra!!! - the current radius
  // appears to be too small
  //  if (enclosingCylinder->ShouldMiss(p, v)) return UUtils::kInfinity;

  // this just takes too much time
  //  idistance = enclosingCylinder->DistanceTo(pb, v);
  //  if (idistance == UUtils::kInfinity) return idistance;

  pb = p + idistance * v;
  int index = GetSection(pb.z());
  pb = p;
  int increment = (v.z() > 0) ? 1 : -1;
  if (std::fabs(v.z()) < fgTolerance) increment = 0;

  double distance = UUtils::kInfinity;
  do
  {
    const UPolyconeSection& section = fSections[index];
    pb.z() -= section.shift;
    distance = section.solid->DistanceToIn(pb, v);
    if (distance < UUtils::kInfinity || !increment)
      break;
    index += increment;
    pb.z() += section.shift;
  }
  while (index >= 0 && index <= fMaxSection);

  return distance;
}


double UPolycone::DistanceToOut(const UVector3&  p, const UVector3& v,
                                UVector3& n, bool& convex, double ) const
{
  UVector3 pn(p);
  
  if (fOriginalParameters->fNumZPlanes == 2)
  {
    const UPolyconeSection& section = fSections[0];
    pn.z() -= section.shift;
    return section.solid->DistanceToOut(pn, v, n, convex);
  }

  int indexLow = GetSection(p.z()-fgTolerance);
  int indexHigh = GetSection(p.z()+fgTolerance);
  int index = 0;
 
  if ( indexLow != indexHigh )  // we are close to Surface,
  {                             // section has to be identified
    const UPolyconeSection& section = fSections[indexLow];
    pn.z() -= section.shift;
    if (section.solid->Inside(pn) == eOutside)
    {
      index=indexHigh;
    }
    else
    {
      index=indexLow;
    }
    pn.z()=p.z();
  }
  else
  {
    index=indexLow;
  } 
  double totalDistance = 0;
  int increment = (v.z() > 0) ? 1 : -1;
  bool convexloc = true;
 
  UVector3 section_norm;
  section_norm.Set(0);
  int istep = 0; 
  do
  {
    const UPolyconeSection& section = fSections[index];
    
    if (totalDistance != 0||(istep < 2))
    {
      pn = p + (totalDistance ) * v; // point must be shifted, so it could
                                     // eventually get into another solid
      pn.z() -= section.shift;
      if (section.solid->Inside(pn) == eOutside)
      {
        break;
      }
    }
    else
    {
      pn.z() -= section.shift;
    }
    istep = istep+1;
    double distance = section.solid-> DistanceToOut(pn, v, n, convexloc);

    if(std::fabs(distance) < 0.5*fgTolerance)  // Section Surface case
    {
      int index1 = index;
      if( ( index > 0) && ( index < fMaxSection ) )
      {
        index1 += increment;
      }
      else
      {
        if( (index == 0) && ( increment > 0 ) )  { index1 += increment; }
        if( (index == fMaxSection) && (increment<0) )  { index1 += increment;}
      }
      UVector3 pte = p+(totalDistance+distance)*v;
      const UPolyconeSection& section1 = fSections[index1];
      pte.z() -= section1.shift;
      if (section1.solid->Inside(pte) == eOutside)
      {
       break;
      }
    }

    if ( (index < fMaxSection) && (index > 0 ) )   // Convexity
    {
      if((convexloc) && (section.convex))
      {
        convexloc = true;
      }
      else
      {
        convexloc = false;
      }
    }
    
    totalDistance += distance;
    index += increment;
  }
  while (index >= 0 && index <= fMaxSection);

  convex=convexloc;
  if (convex)  //Check final convexity for dz
  {
    pn = p + (totalDistance) * v;
    double halfTolerance = 0.5*fgTolerance;
    if( (index < fMaxSection) && (index > 0 ) )
    {
      double dz1 = std::fabs(pn.z()-fZs[index]);
      double dz2 = std::fabs(pn.z()-fZs[index+1]);
      if(dz1 < halfTolerance)  { convex=false; }
      if(dz2 < halfTolerance)  { convex=false; }
    }
    else
    {
      if(index>=fMaxSection)
      {
        if (std::fabs(pn.z()-fZs[fMaxSection]) < halfTolerance) { convex=false; }
      }
      else
      {
        if(index<=0)
        {
          if (std::fabs(pn.z()-fZs[1]) < halfTolerance) { convex=false;}
        }
      }
    }
  }
  convex=false;

  return totalDistance;
}


double UPolycone::SafetyFromInside(const UVector3& p, bool) const
{
  int index = UVoxelizer::BinarySearch(fZs, p.z());
  if (index < 0 || index > fMaxSection) return 0;
  
  double rho=std::sqrt(p.x()*p.x()+p.y()*p.y());
  double safeR = SafetyFromInsideSection(index,rho, p);
  double safeDown = p.z()-fZs[index];
  double safeUp = fZs[index+1]-p.z();
  
  double minSafety =safeR;
  
  if (minSafety == UUtils::kInfinity) return 0;
  if (minSafety < 1e-6) return 0;
 
  for (int i = index + 1; i <= fMaxSection; ++i)
  {
    double dz1 = fZs[i] - p.z();
    double dz2 = fZs[i+1] - p.z();
    safeR = SafetyFromOutsideSection(i,rho, p); 
    if (safeR < 0.) { safeUp=dz1; break; } 
    if (dz1 < dz2) { safeR = std::sqrt(safeR*safeR+dz1*dz1); }
    else {safeR = std::sqrt(safeR*safeR+dz1*dz1); }
    if (safeR < dz1) { safeUp = safeR; break; }
    if (safeR < dz2) { safeUp = safeR; break; }
    safeUp=dz2;
  }

  if (index > 0)
  {
    for (int i = index - 1; i >= 0; --i)
    {
      double dz1 = p.z()-fZs[i+1];
      double dz2 = p.z()-fZs[i];
      safeR = SafetyFromOutsideSection(i,rho, p); 
      if (safeR < 0.) { safeDown=dz1; break; } 
      if(dz1 < dz2) { safeR = std::sqrt(safeR*safeR+dz1*dz1); }
      else { safeR = std::sqrt(safeR*safeR+dz1*dz1); }
      if (safeR < dz1) { safeDown = safeR; break; }
      if (safeR < dz2) { safeDown = safeR; break; }
      safeDown=dz2;
    }
  }
  if (safeUp < minSafety) minSafety=safeUp;
  if (safeDown < minSafety) minSafety=safeDown;
  
  return minSafety;
}

double UPolycone::SafetyFromOutside(const UVector3& p, bool aAccurate) const
{
  if (!aAccurate)
    return enclosingCylinder->SafetyFromOutside(p);

  int index = GetSection(p.z());
  double minSafety = SafetyFromOutsideSection(index, p);
  if (minSafety < 1e-6) return minSafety;

  double zbase = fZs[index + 1];
  for (int i = index + 1; i <= fMaxSection; ++i)
  {
    double dz = fZs[i] - zbase;
    if (dz >= minSafety) break;
    double safety = SafetyFromOutsideSection(i, p);
    if (safety < minSafety) minSafety = safety;
  }

  zbase = fZs[index - 1];
  for (int i = index - 1; i >= 0; --i)
  {
    double dz = zbase - fZs[i];
    if (dz >= minSafety) break;
    double safety = SafetyFromOutsideSection(i, p);
    if (safety < minSafety) minSafety = safety;
  }
  return minSafety;
}

bool UPolycone::Normal(const UVector3& p, UVector3& n) const
{
  double htolerance = 0.5 * fgTolerance;
  int index = GetSection(p.z());

  EnumInside nextPos;
  int nextSection;

  if (index > 0 && p.z()  - fZs[index] < htolerance)
  {
    nextSection = index - 1;
    nextPos = InsideSection(nextSection, p);
  }
  else if (index < fMaxSection && fZs[index + 1] - p.z() < htolerance)
  {
    nextSection = index + 1;
    nextPos = InsideSection(nextSection, p);
  }
  else
  {
    const UPolyconeSection& section = fSections[index];
    UVector3 ps(p.x(), p.y(), p.z() - section.shift);
    bool res = section.solid->Normal(ps, n);
    return res;
  }

  // even if it says we are on the surface, actually it do not have to be

//  "TODO special case when point is on the border of two z-sections",
//    "we should implement this after safety's";

  EnumInside pos = InsideSection(index, p);

  if (nextPos == eInside)
  {
    NormalSection(index, p, n);
    return false;
  }

  if (pos == eSurface && nextPos == eSurface)
  {
    UVector3 n2;
    NormalSection(index, p, n);
    NormalSection(nextSection, p, n2);
    if ((n + n2).Mag2() < 1000 * frTolerance)
    {
      // "we are inside. see TODO above";
      
      NormalSection(index, p, n);
      return false;
    }
  }

  if (nextPos == eSurface || pos == eSurface)
  {
    if (pos != eSurface) index = nextSection;
    bool res = NormalSection(index, p, n);
    return res;
  }

  NormalSection(index, p, n);
  // "we are outside. see TODO above";
  return false;
}

void UPolycone::Extent(UVector3& aMin, UVector3& aMax) const
{
  double r = enclosingCylinder->radius;
  aMin.Set(-r, -r, fZs.front());
  aMax.Set(r, r, fZs.back());
}

double UPolycone::Capacity()
{
  if (fCubicVolume != 0.)
  {
    ;
  }
  else
  {
    for (int i = 0; i < fMaxSection+1; i++)
    {
      UPolyconeSection& section = fSections[i];
      fCubicVolume += section.solid->Capacity();
    }
  }
  return fCubicVolume;
}

double UPolycone::SurfaceArea()
{
  if (fSurfaceArea != 0)
  {
    ;
  }
  else
  {
    double Area = 0, totArea = 0;
    int i = 0;
    int numPlanes = fOriginalParameters->fNumZPlanes;


    std::vector<double> areas;       // (numPlanes+1);
    std::vector<UVector3> points; // (numPlanes-1);

    areas.push_back(UUtils::kPi * (UUtils::sqr(fOriginalParameters->Rmax[0])
                                   - UUtils::sqr(fOriginalParameters->Rmin[0])));

    for (i = 0; i < numPlanes - 1; i++)
    {
      Area = (fOriginalParameters->Rmin[i] + fOriginalParameters->Rmin[i + 1])
             * std::sqrt(UUtils::sqr(fOriginalParameters->Rmin[i]
                                     - fOriginalParameters->Rmin[i + 1]) +
                         UUtils::sqr(fOriginalParameters->fZValues[i + 1]
                                     - fOriginalParameters->fZValues[i]));

      Area += (fOriginalParameters->Rmax[i] + fOriginalParameters->Rmax[i + 1])
              * std::sqrt(UUtils::sqr(fOriginalParameters->Rmax[i]
                                      - fOriginalParameters->Rmax[i + 1]) +
                          UUtils::sqr(fOriginalParameters->fZValues[i + 1]
                                      - fOriginalParameters->fZValues[i]));

      Area *= 0.5 * (endPhi - startPhi);

      if (fOriginalParameters->fOpeningAngle < 2 * UUtils::kPi)
      {
        Area += std::fabs(fOriginalParameters->fZValues[i + 1]
                          - fOriginalParameters->fZValues[i]) *
                (fOriginalParameters->Rmax[i]
                 + fOriginalParameters->Rmax[i + 1]
                 - fOriginalParameters->Rmin[i]
                 - fOriginalParameters->Rmin[i + 1]);
      }
      areas.push_back(Area);
      totArea += Area;
    }

    areas.push_back(UUtils::kPi * (UUtils::sqr(fOriginalParameters->Rmax[numPlanes - 1]) -
                                   UUtils::sqr(fOriginalParameters->Rmin[numPlanes - 1])));

    totArea += (areas[0] + areas[numPlanes]);
    fSurfaceArea = totArea;

  }

  return fSurfaceArea;
}

/////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface
//
// GetPointOnCone
//
// Auxiliary method for Get Point On Surface
//
UVector3 UPolycone::GetPointOnCone(double fRmin1, double fRmax1,
                                   double fRmin2, double fRmax2,
                                   double zOne,   double zTwo,
                                   double& totArea) const
{
  // declare working variables
  //
  double Aone, Atwo, Afive, phi, zRand, fDPhi, cosu, sinu;
  double rRand1, rmin, rmax, chose, rone, rtwo, qone, qtwo;
  double fDz = (zTwo - zOne) / 2., afDz = std::fabs(fDz);
  UVector3 point, offset = UVector3(0., 0., 0.5 * (zTwo + zOne));
  fDPhi = endPhi - startPhi;
  rone = (fRmax1 - fRmax2) / (2.*fDz);
  rtwo = (fRmin1 - fRmin2) / (2.*fDz);
  if (fRmax1 == fRmax2)
  {
    qone = 0.;
  }
  else
  {
    qone = fDz * (fRmax1 + fRmax2) / (fRmax1 - fRmax2);
  }
  if (fRmin1 == fRmin2)
  {
    qtwo = 0.;
  }
  else
  {
    qtwo = fDz * (fRmin1 + fRmin2) / (fRmin1 - fRmin2);
  }
  Aone   = 0.5 * fDPhi * (fRmax2 + fRmax1) * (UUtils::sqr(fRmin1 - fRmin2) + UUtils::sqr(zTwo - zOne));
  Atwo   = 0.5 * fDPhi * (fRmin2 + fRmin1) * (UUtils::sqr(fRmax1 - fRmax2) + UUtils::sqr(zTwo - zOne));
  Afive  = fDz * (fRmax1 - fRmin1 + fRmax2 - fRmin2);
  totArea = Aone + Atwo + 2.*Afive;

  phi  = UUtils::Random(startPhi, endPhi);
  cosu = std::cos(phi);
  sinu = std::sin(phi);


  if ((startPhi == 0) && (endPhi == 2 * UUtils::kPi))
  {
    Afive = 0;
  }
  chose = UUtils::Random(0., Aone + Atwo + 2.*Afive);
  if ((chose >= 0) && (chose < Aone))
  {
    if (fRmax1 != fRmax2)
    {
      zRand = UUtils::Random(-1.*afDz, afDz);
      point = UVector3(rone * cosu * (qone - zRand),
                       rone * sinu * (qone - zRand), zRand);
    }
    else
    {
      point = UVector3(fRmax1 * cosu, fRmax1 * sinu,
                       UUtils::Random(-1.*afDz, afDz));

    }
  }
  else if (chose >= Aone && chose < Aone + Atwo)
  {
    if (fRmin1 != fRmin2)
    {
      zRand = UUtils::Random(-1.*afDz, afDz);
      point = UVector3(rtwo * cosu * (qtwo - zRand),
                       rtwo * sinu * (qtwo - zRand), zRand);

    }
    else
    {
      point = UVector3(fRmin1 * cosu, fRmin1 * sinu,
                       UUtils::Random(-1.*afDz, afDz));
    }
  }
  else if ((chose >= Aone + Atwo + Afive) && (chose < Aone + Atwo + 2.*Afive))
  {
    zRand  = UUtils::Random(-1.*afDz, afDz);
    rmin   = fRmin2 - ((zRand - fDz) / (2.*fDz)) * (fRmin1 - fRmin2);
    rmax   = fRmax2 - ((zRand - fDz) / (2.*fDz)) * (fRmax1 - fRmax2);
    rRand1 = std::sqrt(UUtils::Random() * (UUtils::sqr(rmax) - UUtils::sqr(rmin)) + UUtils::sqr(rmin));
    point  = UVector3(rRand1 * std::cos(startPhi),
                      rRand1 * std::sin(startPhi), zRand);
  }
  else
  {
    zRand  = UUtils::Random(-1.*afDz, afDz);
    rmin   = fRmin2 - ((zRand - fDz) / (2.*fDz)) * (fRmin1 - fRmin2);
    rmax   = fRmax2 - ((zRand - fDz) / (2.*fDz)) * (fRmax1 - fRmax2);
    rRand1 = std::sqrt(UUtils::Random() * (UUtils::sqr(rmax) - UUtils::sqr(rmin)) + UUtils::sqr(rmin));
    point  = UVector3(rRand1 * std::cos(endPhi),
                      rRand1 * std::sin(endPhi), zRand);

  }

  return point + offset;
}


//
// GetPointOnTubs
//
// Auxiliary method for GetPoint On Surface
//
UVector3 UPolycone::GetPointOnTubs(double fRMin, double fRMax,
                                   double zOne,  double zTwo,
                                   double& totArea) const
{
  double xRand, yRand, zRand, phi, cosphi, sinphi, chose,
         aOne, aTwo, aFou, rRand, fDz, fSPhi, fDPhi;
  fDz = std::fabs(0.5 * (zTwo - zOne));
  fSPhi = startPhi;
  fDPhi = endPhi - startPhi;

  aOne = 2.*fDz * fDPhi * fRMax;
  aTwo = 2.*fDz * fDPhi * fRMin;
  aFou = 2.*fDz * (fRMax - fRMin);
  totArea = aOne + aTwo + 2.*aFou;
  phi    = UUtils::Random(startPhi, endPhi);
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);
  rRand  = fRMin + (fRMax - fRMin) * std::sqrt(UUtils::Random());

  if (startPhi == 0 && endPhi == 2 * UUtils::kPi)
    aFou = 0;

  chose  = UUtils::Random(0., aOne + aTwo + 2.*aFou);
  if ((chose >= 0) && (chose < aOne))
  {
    xRand = fRMax * cosphi;
    yRand = fRMax * sinphi;
    zRand = UUtils::Random(-1.*fDz, fDz);
    return UVector3(xRand, yRand, zRand + 0.5 * (zTwo + zOne));
  }
  else if ((chose >= aOne) && (chose < aOne + aTwo))
  {
    xRand = fRMin * cosphi;
    yRand = fRMin * sinphi;
    zRand = UUtils::Random(-1.*fDz, fDz);
    return UVector3(xRand, yRand, zRand + 0.5 * (zTwo + zOne));
  }
  else if ((chose >= aOne + aTwo) && (chose < aOne + aTwo + aFou))
  {
    xRand = rRand * std::cos(fSPhi + fDPhi);
    yRand = rRand * std::sin(fSPhi + fDPhi);
    zRand = UUtils::Random(-1.*fDz, fDz);
    return UVector3(xRand, yRand, zRand + 0.5 * (zTwo + zOne));
  }

  // else

  xRand = rRand * std::cos(fSPhi + fDPhi);
  yRand = rRand * std::sin(fSPhi + fDPhi);
  zRand = UUtils::Random(-1.*fDz, fDz);
  return UVector3(xRand, yRand, zRand + 0.5 * (zTwo + zOne));
}


//
// GetPointOnRing
//
// Auxiliary method for GetPoint On Surface
//
UVector3 UPolycone::GetPointOnRing(double fRMin1, double fRMax1,
                                   double fRMin2, double fRMax2,
                                   double zOne) const
{
  double xRand, yRand, phi, cosphi, sinphi, rRand1, rRand2, A1, Atot, rCh;
  phi    = UUtils::Random(startPhi, endPhi);
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);

  if (fRMin1 == fRMin2)
  {
    rRand1 = fRMin1;
    A1 = 0.;
  }
  else
  {
    rRand1 = UUtils::Random(fRMin1, fRMin2);
    A1 = std::fabs(fRMin2 * fRMin2 - fRMin1 * fRMin1);
  }
  if (fRMax1 == fRMax2)
  {
    rRand2 = fRMax1;
    Atot = A1;
  }
  else
  {
    rRand2 = UUtils::Random(fRMax1, fRMax2);
    Atot   = A1 + std::fabs(fRMax2 * fRMax2 - fRMax1 * fRMax1);
  }
  rCh   = UUtils::Random(0., Atot);

  if (rCh > A1)
  {
    rRand1 = rRand2;
  }

  xRand = rRand1 * cosphi;
  yRand = rRand1 * sinphi;

  return UVector3(xRand, yRand, zOne);
}


//
// GetPointOnCut
//
// Auxiliary method for Get Point On Surface
//
UVector3 UPolycone::GetPointOnCut(double fRMin1, double fRMax1,
                                  double fRMin2, double fRMax2,
                                  double zOne,  double zTwo,
                                  double& totArea) const
{
  if (zOne == zTwo)
  {
    return GetPointOnRing(fRMin1, fRMax1, fRMin2, fRMax2, zOne);
  }
  if ((fRMin1 == fRMin2) && (fRMax1 == fRMax2))
  {
    return GetPointOnTubs(fRMin1, fRMax1, zOne, zTwo, totArea);
  }
  return GetPointOnCone(fRMin1, fRMax1, fRMin2, fRMax2, zOne, zTwo, totArea);
}


//
// GetPointOnSurface
//
UVector3 UPolycone::GetPointOnSurface() const
{
  double Area = 0, totArea = 0, Achose1 = 0, Achose2 = 0, phi, cosphi, sinphi, rRand;
  int i = 0;
  int numPlanes = fOriginalParameters->fNumZPlanes;

  phi = UUtils::Random(startPhi, endPhi);
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);

  rRand = fOriginalParameters->Rmin[0] +
          ((fOriginalParameters->Rmax[0] - fOriginalParameters->Rmin[0])
           * std::sqrt(UUtils::Random()));

  std::vector<double> areas;       // (numPlanes+1);
  std::vector<UVector3> points; // (numPlanes-1);

  areas.push_back(UUtils::kPi * (UUtils::sqr(fOriginalParameters->Rmax[0])
                                 - UUtils::sqr(fOriginalParameters->Rmin[0])));

  for (i = 0; i < numPlanes - 1; i++)
  {
    Area = (fOriginalParameters->Rmin[i] + fOriginalParameters->Rmin[i + 1])
           * std::sqrt(UUtils::sqr(fOriginalParameters->Rmin[i]
                                   - fOriginalParameters->Rmin[i + 1]) +
                       UUtils::sqr(fOriginalParameters->fZValues[i + 1]
                                   - fOriginalParameters->fZValues[i]));

    Area += (fOriginalParameters->Rmax[i] + fOriginalParameters->Rmax[i + 1])
            * std::sqrt(UUtils::sqr(fOriginalParameters->Rmax[i]
                                    - fOriginalParameters->Rmax[i + 1]) +
                        UUtils::sqr(fOriginalParameters->fZValues[i + 1]
                                    - fOriginalParameters->fZValues[i]));

    Area *= 0.5 * (endPhi - startPhi);

    if (startPhi == 0. && endPhi == 2 * UUtils::kPi)
    {
      Area += std::fabs(fOriginalParameters->fZValues[i + 1]
                        - fOriginalParameters->fZValues[i]) *
              (fOriginalParameters->Rmax[i]
               + fOriginalParameters->Rmax[i + 1]
               - fOriginalParameters->Rmin[i]
               - fOriginalParameters->Rmin[i + 1]);
    }
    areas.push_back(Area);
    totArea += Area;
  }

  areas.push_back(UUtils::kPi * (UUtils::sqr(fOriginalParameters->Rmax[numPlanes - 1]) -
                                 UUtils::sqr(fOriginalParameters->Rmin[numPlanes - 1])));

  totArea += (areas[0] + areas[numPlanes]);
  double chose = UUtils::Random(0., totArea);

  if ((chose >= 0.) && (chose < areas[0]))
  {
    return UVector3(rRand * cosphi, rRand * sinphi,
                    fOriginalParameters->fZValues[0]);
  }

  for (i = 0; i < numPlanes - 1; i++)
  {
    Achose1 += areas[i];
    Achose2 = (Achose1 + areas[i + 1]);
    if (chose >= Achose1 && chose < Achose2)
    {
      return GetPointOnCut(fOriginalParameters->Rmin[i],
                           fOriginalParameters->Rmax[i],
                           fOriginalParameters->Rmin[i + 1],
                           fOriginalParameters->Rmax[i + 1],
                           fOriginalParameters->fZValues[i],
                           fOriginalParameters->fZValues[i + 1], Area);
    }
  }

  rRand = fOriginalParameters->Rmin[numPlanes - 1] +
          ((fOriginalParameters->Rmax[numPlanes - 1] - fOriginalParameters->Rmin[numPlanes - 1])
           * std::sqrt(UUtils::Random()));

  return UVector3(rRand * cosphi, rRand * sinphi,
                  fOriginalParameters->fZValues[numPlanes - 1]);

}

//
// UPolyconeHistorical stuff
//

UPolyconeHistorical::UPolyconeHistorical()
  : fStartAngle(0.), fOpeningAngle(0.), fNumZPlanes(0),
    fZValues(0), Rmin(0), Rmax(0)
{
}

UPolyconeHistorical::~UPolyconeHistorical()
{
}

UPolyconeHistorical::
UPolyconeHistorical(const UPolyconeHistorical& source)
{
  fStartAngle  = source.fStartAngle;
  fOpeningAngle = source.fOpeningAngle;
  fNumZPlanes = source.fNumZPlanes;

  fZValues.resize(fNumZPlanes);
  Rmin.resize(fNumZPlanes);
  Rmax.resize(fNumZPlanes);

  for (int i = 0; i < fNumZPlanes; i++)
  {
    fZValues[i] = source.fZValues[i];
    Rmin[i]    = source.Rmin[i];
    Rmax[i]    = source.Rmax[i];
  }
}

UPolyconeHistorical&
UPolyconeHistorical::operator=(const UPolyconeHistorical& right)
{
  if (&right == this) return *this;

  fStartAngle  = right.fStartAngle;
  fOpeningAngle = right.fOpeningAngle;
  fNumZPlanes = right.fNumZPlanes;

  fZValues.resize(fNumZPlanes);
  Rmin.resize(fNumZPlanes);
  Rmax.resize(fNumZPlanes);

  for (int i = 0; i < fNumZPlanes; i++)
  {
    fZValues[i] = right.fZValues[i];
    Rmin[i]    = right.Rmin[i];
    Rmax[i]    = right.Rmax[i];
  }

  return *this;
}

VUSolid* UPolycone::Clone() const
{
  return new UPolycone(*this);
}
//
// Copy constructor
//
UPolycone::UPolycone(const UPolycone& source): VUSolid(source)
{
  CopyStuff(source);
}


//
// Assignment operator
//
UPolycone& UPolycone::operator=(const UPolycone& source)
{
  if (this == &source) return *this;

  //VUSolid::operator=( source );

  //delete [] corners;

  delete enclosingCylinder;

  CopyStuff(source);

  return *this;
}
//
// CopyStuff
//
void UPolycone::CopyStuff(const UPolycone& source)
{
  //
  // Simple stuff
  //

  startPhi  = source.startPhi;
  endPhi    = source.endPhi;
  phiIsOpen = source.phiIsOpen;
  fCubicVolume    = source.fCubicVolume;
  fSurfaceArea    = source.fSurfaceArea;
  fBox = source.fBox;
  //
  // The array of planes
  //
  fOriginalParameters = source.fOriginalParameters;
  //
  // Enclosing cylinder
  //
  enclosingCylinder = new UEnclosingCylinder(*source.enclosingCylinder);
}
//
// Get Entity Type
//
UGeometryType UPolycone::GetEntityType() const
{
      return "Polycone";
}
//
// Set Original Parameters
//
bool  UPolycone::SetOriginalParameters(UReduciblePolygon* rz)
{
  int numPlanes = (int)numCorner;
  bool isConvertible = true;
  double Zmax = rz->Bmax();
  rz->StartWithZMin();

  // Prepare vectors for storage
  //
  std::vector<double> Z;
  std::vector<double> Rmin;
  std::vector<double> Rmax;

  int countPlanes = 1;
  int icurr = 0;
  int icurl = 0;

  // first plane Z=Z[0]
  //
  Z.push_back(corners[0].z);
  double Zprev = Z[0];
  if (Zprev == corners[1].z)
  {
    Rmin.push_back(corners[0].r);
    Rmax.push_back(corners[1].r);
    icurr = 1;
  }
  else if (Zprev == corners[numPlanes - 1].z)
  {
    Rmin.push_back(corners[numPlanes - 1].r);
    Rmax.push_back(corners[0].r);
    icurl = numPlanes - 1;
  }
  else
  {
    Rmin.push_back(corners[0].r);
    Rmax.push_back(corners[0].r);
  }

  // next planes until last
  //
  int inextr = 0, inextl = 0;
  for (int i = 0; i < numPlanes - 2; i++)
  {
    inextr = 1 + icurr;
    inextl = (icurl <= 0) ? numPlanes - 1 : icurl - 1;

    if ((corners[inextr].z >= Zmax) & (corners[inextl].z >= Zmax))
    {
      break;
    }

    double Zleft = corners[inextl].z;
    double Zright = corners[inextr].z;
    if (Zright > Zleft)
    {
      Z.push_back(Zleft);
      countPlanes++;
      double difZr = corners[inextr].z - corners[icurr].z;
      double difZl = corners[inextl].z - corners[icurl].z;

      if (std::fabs(difZl) < frTolerance)
      {
        if (std::fabs(difZr) < frTolerance)
        {
          Rmin.push_back(corners[inextl].r);
          Rmax.push_back(corners[icurr].r);
        }
        else
        {
          Rmin.push_back(corners[inextl].r);
          Rmax.push_back(corners[icurr].r + (Zleft-corners[icurr].z)/difZr
                                *(corners[inextr].r - corners[icurr].r)); 
        }
      }
      else if (difZl >= frTolerance)
      {
        if (std::fabs(difZr) < frTolerance)
        {
          Rmin.push_back(corners[icurl].r);
          Rmax.push_back(corners[icurr].r);
        }
        else
        {
          Rmin.push_back(corners[icurl].r);
          Rmax.push_back(corners[icurr].r + (Zleft-corners[icurr].z)/difZr
                                *(corners[inextr].r - corners[icurr].r));
        }
      }
      else
      {
        isConvertible = false;
        break;
      }
      icurl = (icurl == 0) ? numPlanes - 1 : icurl - 1;
    }
    else if (std::fabs(Zright - Zleft) < frTolerance) // Zright=Zleft
    {
      Z.push_back(Zleft);
      countPlanes++;
      icurr++;

      icurl = (icurl == 0) ? numPlanes - 1 : icurl - 1;

      Rmin.push_back(corners[inextl].r);
      Rmax.push_back(corners[inextr].r);
    }
    else  // Zright<Zleft
    {
      Z.push_back(Zright);
      countPlanes++;

      double difZr = corners[inextr].z - corners[icurr].z;
      double difZl = corners[inextl].z - corners[icurl].z;
      if (std::fabs(difZr) < frTolerance)
      {
        if (std::fabs(difZl) < frTolerance)
        {
          Rmax.push_back(corners[inextr].r);
          Rmin.push_back (corners[icurr].r); 
        } 
        else
        {
          Rmin.push_back (corners[icurl].r + (Zright-corners[icurl].z)/difZl
                                 * (corners[inextl].r - corners[icurl].r));
          Rmax.push_back(corners[inextr].r);
        }
        icurr++;
      }           // plate
      else if (difZr >= frTolerance)
      {
        if (std::fabs(difZl) < frTolerance)
        {
          Rmax.push_back(corners[inextr].r);
          Rmin.push_back(corners[icurr].r);
        }
        else
        {
          Rmax.push_back(corners[inextr].r);
          Rmin.push_back(corners[icurl].r + (Zright - corners[icurl].z) / difZl
                         * (corners[inextl].r - corners[icurl].r));
        }
        icurr++;
      }
      else
      {
        isConvertible = false;
        break;
      }
    }
  }   // end for loop

  // last plane Z=Zmax
  //
  Z.push_back(Zmax);
  countPlanes++;
  inextr = 1 + icurr;
  inextl = (icurl <= 0) ? numPlanes - 1 : icurl - 1;

  if (corners[inextr].z == corners[inextl].z)
  {
    Rmax.push_back(corners[inextr].r);
    Rmin.push_back(corners[inextl].r);
  }
  else
  {
    Rmax.push_back(corners[inextr].r);
    Rmin.push_back(corners[inextl].r);
  }

  // Set original parameters Rmin,Rmax,Z
  //
  if (isConvertible)
  {
    fOriginalParameters = new UPolyconeHistorical;
    fOriginalParameters->fZValues.resize(numPlanes);
    fOriginalParameters->Rmin.resize(numPlanes);
    fOriginalParameters->Rmax.resize(numPlanes);

    for (int j = 0; j < countPlanes; j++)
    {
      fOriginalParameters->fZValues[j] = Z[j];
      fOriginalParameters->Rmax[j] = Rmax[j];
      fOriginalParameters->Rmin[j] = Rmin[j];
    }
    fOriginalParameters->fStartAngle = startPhi;
    fOriginalParameters->fOpeningAngle = endPhi - startPhi;
    fOriginalParameters->fNumZPlanes = countPlanes;

  }
  else  // Set parameters(r,z) with Rmin==0 as convention
  {
    std::ostringstream message;
    message << "Polycone " << GetName() << std::endl
            << "cannot be converted to Polycone with (Rmin,Rmaz,Z) parameters!";
    UUtils::Exception("UPolycone::SetOriginalParameters()", "GeomSolids0002",
                      UWarning, 1, "can not convert");

    fOriginalParameters = new UPolyconeHistorical;

    fOriginalParameters->fZValues.resize(numPlanes);
    fOriginalParameters->Rmin.resize(numPlanes);
    fOriginalParameters->Rmax.resize(numPlanes);

    for (int j = 0; j < numPlanes; j++)
    {
      fOriginalParameters->fZValues[j] = corners[j].z;
      fOriginalParameters->Rmax[j] = corners[j].r;
      fOriginalParameters->Rmin[j] = 0.0;
    }
    fOriginalParameters->fStartAngle = startPhi;
    fOriginalParameters->fOpeningAngle = endPhi - startPhi;
    fOriginalParameters->fNumZPlanes = numPlanes;
  }
  return isConvertible;
}
//
// Reset 
//
void UPolycone::Reset()
{
   //
   // Clear old setup
   //
   delete enclosingCylinder;
 
  fCubicVolume = 0;
  fSurfaceArea = 0;
  double phiStart=fOriginalParameters->fStartAngle;
  double* Z, *R1, *R2;
  int num = fOriginalParameters->fNumZPlanes;
    Z = new double[num];
    R1 = new double[num];
    R2 = new double[num];
    for (int i = 0; i < num; i++)
    {
      Z[i] = fOriginalParameters->fZValues[i];
      R1[i] = fOriginalParameters->Rmin[i];
      R2[i] = fOriginalParameters->Rmax[i];
    }

   Init(phiStart, phiStart+ fOriginalParameters->fOpeningAngle, num, Z, R1, R2);
    delete [] R1;
    delete [] Z;
    delete [] R2;
}
