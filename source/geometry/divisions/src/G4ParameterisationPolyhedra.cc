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
// $Id: G4ParameterisationPolyhedra.cc 92635 2015-09-09 13:39:54Z gcosmo $
//
// class G4ParameterisationPolyhedra Implementation file
//
// 14.10.03 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// --------------------------------------------------------------------

#include "G4ParameterisationPolyhedra.hh"

#include <iomanip>
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4GeometryTolerance.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ReflectedSolid.hh"
#include "G4Polyhedra.hh"

//--------------------------------------------------------------------------
G4VParameterisationPolyhedra::
G4VParameterisationPolyhedra( EAxis axis, G4int nDiv, G4double width,
                              G4double offset, G4VSolid* msolid,
                              DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
  std::ostringstream message;
#ifdef G4MULTITHREADED
  message << "Divisions for G4Polyhedra currently NOT supported in MT-mode."
          << G4endl
          << "Sorry! Solid: " << msolid->GetName();
  G4Exception("G4VParameterisationPolyhedra::G4VParameterisationPolyhedra()",
              "GeomDiv0001", FatalException, message);
#endif

  G4Polyhedra* msol = (G4Polyhedra*)(msolid);
  if ((msolid->GetEntityType() != "G4ReflectedSolid") && (msol->IsGeneric()))
  {
    message << "Generic construct for G4Polyhedra NOT supported." << G4endl
            << "Sorry! Solid: " << msol->GetName();
    G4Exception("G4VParameterisationPolyhedra::G4VParameterisationPolyhedra()",
                "GeomDiv0001", FatalException, message);
  }
  if (msolid->GetEntityType() == "G4ReflectedSolid")
  {
    // Get constituent solid  
    G4VSolid* mConstituentSolid 
       = ((G4ReflectedSolid*)msolid)->GetConstituentMovedSolid();
    msol = (G4Polyhedra*)(mConstituentSolid);
  
    // Get parameters
    G4int   nofSides = msol->GetOriginalParameters()->numSide;
    G4int   nofZplanes = msol->GetOriginalParameters()->Num_z_planes;
    G4double* zValues  = msol->GetOriginalParameters()->Z_values;
    G4double* rminValues  = msol->GetOriginalParameters()->Rmin;
    G4double* rmaxValues  = msol->GetOriginalParameters()->Rmax;

    // Invert z values, 
    // convert radius parameters
    G4double* rminValues2 = new G4double[nofZplanes];
    G4double* rmaxValues2 = new G4double[nofZplanes];
    G4double* zValuesRefl = new G4double[nofZplanes];
    for (G4int i=0; i<nofZplanes; i++)
    {
      rminValues2[i] = rminValues[i] * ConvertRadiusFactor(*msol);
      rmaxValues2[i] = rmaxValues[i] * ConvertRadiusFactor(*msol);
      zValuesRefl[i] = - zValues[i];
    }  
    
    G4Polyhedra* newSolid
      = new G4Polyhedra(msol->GetName(),
                        msol->GetStartPhi(), 
                        msol->GetEndPhi() - msol->GetStartPhi(),
                        nofSides,
                        nofZplanes, zValuesRefl, rminValues2, rmaxValues2);

    delete [] rminValues2;       
    delete [] rmaxValues2;       
    delete [] zValuesRefl;       

    msol = newSolid;
    fmotherSolid = newSolid;
    fReflectedSolid = true;
    fDeleteSolid = true;
  }    
}

//------------------------------------------------------------------------
G4VParameterisationPolyhedra::~G4VParameterisationPolyhedra()
{
}

//--------------------------------------------------------------------------
G4double 
G4VParameterisationPolyhedra::
ConvertRadiusFactor(const G4Polyhedra& phedra) const
{
  G4double phiTotal = phedra.GetEndPhi() - phedra.GetStartPhi();
  G4int nofSides = phedra.GetOriginalParameters()->numSide;
  
  if ( (phiTotal <=0) || (phiTotal >
        2*pi+G4GeometryTolerance::GetInstance()->GetAngularTolerance()) )
    { phiTotal = 2*pi; }
  
  return std::cos(0.5*phiTotal/nofSides);
}  

//--------------------------------------------------------------------------
G4ParameterisationPolyhedraRho::
G4ParameterisationPolyhedraRho( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationPolyhedra( axis, nDiv, width, offset, msolid, divType )
{

  CheckParametersValidity();
  SetType( "DivisionPolyhedraRho" );

  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);
  G4PolyhedraHistorical* original_pars = msol->GetOriginalParameters();

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( original_pars->Rmax[0]
                         - original_pars->Rmin[0], width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( original_pars->Rmax[0]
                           - original_pars->Rmin[0], nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationPolyhedraRho - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationPolyhedraRho::~G4ParameterisationPolyhedraRho()
{
}

//---------------------------------------------------------------------
void G4ParameterisationPolyhedraRho::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();

  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);

  if( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH )
  {
    std::ostringstream message;
    message << "In solid " << msol->GetName() << G4endl
            << "Division along R will be done with a width "
            << "different for each solid section." << G4endl
            << "WIDTH will not be used !";
    G4Exception("G4ParameterisationPolyhedraRho::CheckParametersValidity()",
                "GeomDiv1001", JustWarning, message);
  }
  if( foffset != 0. )
  {
    std::ostringstream message;
    message << "In solid " << msol->GetName() << G4endl
            << "Division along  R will be done with a width "
            << "different for each solid section." << G4endl
            << "OFFSET will not be used !";
    G4Exception("G4ParameterisationPolyhedraRho::CheckParametersValidity()",
                "GeomDiv1001", JustWarning, message);
  }
}

//------------------------------------------------------------------------
G4double G4ParameterisationPolyhedraRho::GetMaxParameter() const
{
  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);
  G4PolyhedraHistorical* original_pars = msol->GetOriginalParameters();
  return original_pars->Rmax[0] - original_pars->Rmin[0];
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyhedraRho::
ComputeTransformation( const G4int, G4VPhysicalVolume* physVol ) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.);

  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix: unit

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyhedraRho " << G4endl
           << " foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
#endif

  ChangeRotMatrix( physVol );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyhedraRho "
           << G4endl
           << " Position: " << origin
           << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
#endif
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyhedraRho::
ComputeDimensions( G4Polyhedra& phedra, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);

  G4PolyhedraHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyhedraHistorical origparam( *origparamMother );
  G4int nZplanes = origparamMother->Num_z_planes;

  G4double width = 0.;
  for( G4int ii = 0; ii < nZplanes; ii++ )
  {
    width = CalculateWidth( origparamMother->Rmax[ii]
                          - origparamMother->Rmin[ii], fnDiv, foffset );
    origparam.Rmin[ii] = origparamMother->Rmin[ii]+foffset+width*copyNo;
    origparam.Rmax[ii] = origparamMother->Rmin[ii]+foffset+width*(copyNo+1);
  }

  phedra.SetOriginalParameters(&origparam); // copy values & transfer pointers
  phedra.Reset();                           // reset to new solid parameters

#ifdef G4DIVDEBUG
  if( verbose >= -2 )
  {
    G4cout << "G4ParameterisationPolyhedraRho::ComputeDimensions()" << G4endl
           << "-- Parametrised phedra copy-number: " << copyNo << G4endl;
    phedra.DumpInfo();
  }
#endif
}

//--------------------------------------------------------------------------
G4ParameterisationPolyhedraPhi::
G4ParameterisationPolyhedraPhi( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationPolyhedra( axis, nDiv, width, offset, msolid, divType )
{ 
  CheckParametersValidity();
  SetType( "DivisionPolyhedraPhi" );

  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);
  G4double deltaPhi = msol->GetEndPhi() - msol->GetStartPhi();

  if( divType == DivWIDTH )
  {
    fnDiv = msol->GetNumSide();
  }

  fwidth = CalculateWidth( deltaPhi, fnDiv, 0.0 );

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationPolyhedraPhi - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//------------------------------------------------------------------------
G4ParameterisationPolyhedraPhi::~G4ParameterisationPolyhedraPhi()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationPolyhedraPhi::GetMaxParameter() const
{
  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);
  return msol->GetEndPhi() - msol->GetStartPhi();
}

//---------------------------------------------------------------------
void G4ParameterisationPolyhedraPhi::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();

  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);

  if( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH )
  {
    std::ostringstream message;
    message << "In solid " << msol->GetName() << G4endl
            << " Division along PHI will be done splitting "
            << "in the defined numSide." << G4endl
            << "WIDTH will not be used !";
    G4Exception("G4ParameterisationPolyhedraPhi::CheckParametersValidity()",
                "GeomDiv1001", JustWarning, message);
  }
  if( foffset != 0. )
  {
    std::ostringstream message;
    message << "In solid " << msol->GetName() << G4endl
            << "Division along PHI will be done splitting "
            << "in the defined numSide." << G4endl
            << "OFFSET will not be used !";
    G4Exception("G4ParameterisationPolyhedraPhi::CheckParametersValidity()",
                "GeomDiv1001", JustWarning, message);
  }

  G4PolyhedraHistorical* origparamMother = msol->GetOriginalParameters();

  if( origparamMother->numSide != fnDiv &&  fDivisionType != DivWIDTH)
  { 
    std::ostringstream message;
    message << "Configuration not supported." << G4endl
            << "Division along PHI will be done splitting in the defined"
            << G4endl
            << "numSide, i.e, the number of division would be :"
            << origparamMother->numSide << " instead of " << fnDiv << " !"; 
    G4Exception("G4ParameterisationPolyhedraPhi::CheckParametersValidity()",
                "GeomDiv0001", FatalException, message);
  }
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyhedraPhi::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume *physVol ) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix (so that all volumes point to the centre)
  G4double posi = copyNo*fwidth;

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyhedraPhi - position: " << posi/deg
           << G4endl
           << " copyNo: " << copyNo
           << " - fwidth: " << fwidth/deg << G4endl;
  }
#endif

  ChangeRotMatrix( physVol, -posi );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyhedraPhi " << copyNo
           << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
#endif
}

//--------------------------------------------------------------------------
void
G4ParameterisationPolyhedraPhi::
ComputeDimensions( G4Polyhedra& phedra, const G4int,
                   const G4VPhysicalVolume* ) const
{
  G4Polyhedra* msol = (G4Polyhedra*)(fmotherSolid);

  G4PolyhedraHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyhedraHistorical origparam( *origparamMother );

  origparam.numSide = 1;
  origparam.Start_angle = origparamMother->Start_angle;
  origparam.Opening_angle = fwidth;

  phedra.SetOriginalParameters(&origparam);  // copy values & transfer pointers
  phedra.Reset();                            // reset to new solid parameters

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << "G4ParameterisationPolyhedraPhi::ComputeDimensions():" << G4endl;
    phedra.DumpInfo();
  }
#endif
}

//--------------------------------------------------------------------------
G4ParameterisationPolyhedraZ::
G4ParameterisationPolyhedraZ( EAxis axis, G4int nDiv,
                             G4double width, G4double offset,
                             G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationPolyhedra( axis, nDiv, width, offset, msolid, divType ),
     fNSegment(0),
     fOrigParamMother(((G4Polyhedra*)fmotherSolid)->GetOriginalParameters())
{ 
  CheckParametersValidity();
  SetType( "DivisionPolyhedraZ" );

  if( divType == DivWIDTH )
    {
    fnDiv =
      CalculateNDiv( fOrigParamMother->Z_values[fOrigParamMother->Num_z_planes-1]
                     - fOrigParamMother->Z_values[0] , width, offset );
  }
  else if( divType == DivNDIV )
    {
    fwidth =
      CalculateNDiv( fOrigParamMother->Z_values[fOrigParamMother->Num_z_planes-1]
                     - fOrigParamMother->Z_values[0] , nDiv, offset );
  }
  
#ifdef G4DIVDEBUG
  if( verbose >= 1 )
    {
    G4cout << " G4ParameterisationPolyhedraZ - # divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//---------------------------------------------------------------------
G4ParameterisationPolyhedraZ::~G4ParameterisationPolyhedraZ()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationPolyhedraZ::GetR(G4double z, 
                                           G4double z1, G4double r1,  
                                           G4double z2, G4double r2) const
{
  // Linear parameterisation: 
  // r = az + b
  // a = (r1 - r2)/(z1-z2)
  // b = r1 - a*z1

  return (r1-r2)/(z1-z2)*z + ( r1 - (r1-r2)/(z1-z2)*z1 ) ;
}  
                                           
//------------------------------------------------------------------------
G4double G4ParameterisationPolyhedraZ::GetRmin(G4double z, G4int nseg) const
{
// Get Rmin in the given z position for the given polyhedra segment 

  return GetR(z, 
              fOrigParamMother->Z_values[nseg], 
              fOrigParamMother->Rmin[nseg],
              fOrigParamMother->Z_values[nseg+1], 
              fOrigParamMother->Rmin[nseg+1]);
}  
                                           
//------------------------------------------------------------------------
G4double G4ParameterisationPolyhedraZ::GetRmax(G4double z, G4int nseg) const
{
// Get Rmax in the given z position for the given polyhedra segment 

  return GetR(z, 
              fOrigParamMother->Z_values[nseg], 
              fOrigParamMother->Rmax[nseg],
              fOrigParamMother->Z_values[nseg+1], 
              fOrigParamMother->Rmax[nseg+1]);
}  
                                           
//------------------------------------------------------------------------
G4double G4ParameterisationPolyhedraZ::GetMaxParameter() const
{
  return std::abs (fOrigParamMother->Z_values[fOrigParamMother->Num_z_planes-1]
             -fOrigParamMother->Z_values[0]);
}

//---------------------------------------------------------------------
void G4ParameterisationPolyhedraZ::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();

  // Division will be following the mother polyhedra segments
  if( fDivisionType == DivNDIV ) {
    if( fOrigParamMother->Num_z_planes-1 != fnDiv ) { 
      std::ostringstream message;
      message << "Configuration not supported." << G4endl
              << "Division along Z will be done splitting in the defined"
              << G4endl
              << "Z planes, i.e, the number of division would be :"
              << fOrigParamMother->Num_z_planes-1 << " instead of "
              << fnDiv << " !"; 
      G4Exception("G4ParameterisationPolyhedraZ::CheckParametersValidity()",
                  "GeomDiv0001", FatalException, message);
    }
  }  

  // Division will be done within one polyhedra segment
  // with applying given width and offset
  if( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH ) {
    // Check if divided region does not span over more
    // than one z segment

    G4int isegstart = -1;  // number of the segment containing start position
    G4int isegend = -1;    // number of the segment containing end position

    if ( ! fReflectedSolid ) {
      // The start/end position of the divided region
      G4double zstart 
        = fOrigParamMother->Z_values[0] + foffset;
      G4double zend 
        = fOrigParamMother->Z_values[0] + foffset + fnDiv* fwidth;
   
      G4int counter = 0;
      while ( isegend < 0 && counter < fOrigParamMother->Num_z_planes-1 ) {
        // first segment
        if ( zstart >= fOrigParamMother->Z_values[counter]  &&
             zstart  < fOrigParamMother->Z_values[counter+1] ) {
           isegstart = counter;
        }     
        // last segment
        if ( zend  > fOrigParamMother->Z_values[counter] &&
             zend <= fOrigParamMother->Z_values[counter+1] ) {
           isegend = counter;
        }   
        ++counter;   
      }  // Loop checking, 06.08.2015, G.Cosmo
    }
    else  {
      // The start/end position of the divided region
      G4double zstart 
        = fOrigParamMother->Z_values[0] - foffset;
      G4double zend 
        = fOrigParamMother->Z_values[0] - ( foffset + fnDiv* fwidth);
   
      G4int counter = 0;
      while ( isegend < 0 && counter < fOrigParamMother->Num_z_planes-1 ) {
        // first segment
        if ( zstart <= fOrigParamMother->Z_values[counter]  &&
             zstart  > fOrigParamMother->Z_values[counter+1] ) {
           isegstart = counter;
        }     
        // last segment
        if ( zend  < fOrigParamMother->Z_values[counter] &&
             zend >= fOrigParamMother->Z_values[counter+1] ) {
           isegend = counter;
        }   
        ++counter;   
      }  // Loop checking, 06.08.2015, G.Cosmo
    }
  
    if ( isegstart != isegend ) {
      std::ostringstream message;
      message << "Configuration not supported." << G4endl
              << "Division with user defined width." << G4endl
              << "Solid " << fmotherSolid->GetName() << G4endl
              << "Divided region is not between two Z planes."; 
      G4Exception("G4ParameterisationPolyhedraZ::CheckParametersValidity()",
                  "GeomDiv0001", FatalException, message);
    }
  
    fNSegment = isegstart;
  }  
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyhedraZ::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  if ( fDivisionType == DivNDIV ) {
    // The position of the centre of copyNo-th mother polycone segment
    G4double posi = ( fOrigParamMother->Z_values[copyNo]
                    + fOrigParamMother->Z_values[copyNo+1])/2;
    physVol->SetTranslation( G4ThreeVector(0, 0, posi) );
  }
  
  if ( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH ) {
    // The position of the centre of copyNo-th division

    G4double posi = fOrigParamMother->Z_values[0];
    
    if ( ! fReflectedSolid )
      posi += foffset + (2*copyNo + 1) * fwidth/2.;
    else
      posi -= foffset + (2*copyNo + 1) * fwidth/2.;
    
    physVol->SetTranslation( G4ThreeVector(0, 0, posi) );
  }   

  //----- calculate rotation matrix: unit

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyhedraZ - position: " << posi << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
#endif

  ChangeRotMatrix( physVol );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyhedraZ "
           << copyNo << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
#endif
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyhedraZ::
ComputeDimensions( G4Polyhedra& phedra, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  // Define division solid
  G4PolyhedraHistorical origparam;
  G4int nz = 2; 
  origparam.Num_z_planes = nz;
  origparam.numSide = fOrigParamMother->numSide;
  origparam.Start_angle = fOrigParamMother->Start_angle;
  origparam.Opening_angle = fOrigParamMother->Opening_angle;

  // Define division solid z sections
  origparam.Z_values = new G4double[nz];
  origparam.Rmin = new G4double[nz];
  origparam.Rmax = new G4double[nz];
  origparam.Z_values[0] = - fwidth/2.;
  origparam.Z_values[1] = fwidth/2.;

  if ( fDivisionType == DivNDIV ) {
    // The position of the centre of copyNo-th mother polycone segment
    G4double posi = ( fOrigParamMother->Z_values[copyNo]
                    + fOrigParamMother->Z_values[copyNo+1])/2;

    origparam.Z_values[0] = fOrigParamMother->Z_values[copyNo] - posi;
    origparam.Z_values[1] = fOrigParamMother->Z_values[copyNo+1] - posi;
    origparam.Rmin[0] = fOrigParamMother->Rmin[copyNo];
    origparam.Rmin[1] = fOrigParamMother->Rmin[copyNo+1];
    origparam.Rmax[0] = fOrigParamMother->Rmax[copyNo];
    origparam.Rmax[1] = fOrigParamMother->Rmax[copyNo+1];
  }  

  if ( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH ) {
    if ( ! fReflectedSolid ) {
      origparam.Z_values[0] = - fwidth/2.;
      origparam.Z_values[1] = fwidth/2.;

      // The position of the centre of copyNo-th division
      G4double posi 
        = fOrigParamMother->Z_values[0] + foffset + (2*copyNo + 1) * fwidth/2.;
    
      // The first and last z sides z values
      G4double zstart = posi - fwidth/2.;
      G4double zend = posi + fwidth/2.;
      origparam.Rmin[0] = GetRmin(zstart, fNSegment); 
      origparam.Rmax[0] = GetRmax(zstart, fNSegment);  
      origparam.Rmin[1] = GetRmin(zend, fNSegment); 
      origparam.Rmax[1] = GetRmax(zend, fNSegment); 
    }
    else {   
      origparam.Z_values[0] = fwidth/2.;
      origparam.Z_values[1] = - fwidth/2.;

      // The position of the centre of copyNo-th division
      G4double posi 
        = fOrigParamMother->Z_values[0] - ( foffset + (2*copyNo + 1) * fwidth/2.);
    
      // The first and last z sides z values
      G4double zstart = posi + fwidth/2.;
      G4double zend = posi - fwidth/2.;
      origparam.Rmin[0] = GetRmin(zstart, fNSegment); 
      origparam.Rmax[0] = GetRmax(zstart, fNSegment);  
      origparam.Rmin[1] = GetRmin(zend, fNSegment); 
      origparam.Rmax[1] = GetRmax(zend, fNSegment); 
    } 

    // It can happen due to rounding errors
    if ( origparam.Rmin[0]    < 0.0 ) origparam.Rmin[0] = 0.0;
    if ( origparam.Rmin[nz-1] < 0.0 ) origparam.Rmin[1] = 0.0;
  }  

  phedra.SetOriginalParameters(&origparam);  // copy values & transfer pointers
  phedra.Reset();                            // reset to new solid parameters

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << "G4ParameterisationPolyhedraZ::ComputeDimensions()" << G4endl
           << "-- Parametrised phedra copy-number: " << copyNo << G4endl;
    phedra.DumpInfo();
  }
#endif
}
