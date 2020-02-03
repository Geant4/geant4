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
// G4ParameterisationPolycone[Rho/Phi/Z] implementation
//
// 26.05.03 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
//---------------------------------------------------------------------

#include "G4ParameterisationPolycone.hh"

#include <iomanip>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ReflectedSolid.hh"

//-----------------------------------------------------------------------
G4VParameterisationPolycone::
G4VParameterisationPolycone( EAxis axis, G4int nDiv, G4double width,
                             G4double offset, G4VSolid* msolid,
                             DivisionType divType )
  :  G4VDivisionParameterisation( axis, nDiv, width, offset, divType, msolid )
{
#ifdef G4MULTITHREADED
   std::ostringstream message;
   message << "Divisions for G4Polycone currently NOT supported in MT-mode."
           << G4endl
           << "Sorry! Solid: " << msolid->GetName();
   G4Exception("G4VParameterisationPolycone::G4VParameterisationPolycone()",
               "GeomDiv0001", FatalException, message);
#endif
  G4Polycone* msol = (G4Polycone*)(msolid);
  if (msolid->GetEntityType() == "G4ReflectedSolid")
  {
    // Get constituent solid
    //
    G4VSolid* mConstituentSolid 
       = ((G4ReflectedSolid*)msolid)->GetConstituentMovedSolid();
    msol = (G4Polycone*)(mConstituentSolid);
  
    // Get parameters
    //
    G4int   nofZplanes = msol->GetOriginalParameters()->Num_z_planes;
    G4double* zValues  = msol->GetOriginalParameters()->Z_values;
    G4double* rminValues  = msol->GetOriginalParameters()->Rmin;
    G4double* rmaxValues  = msol->GetOriginalParameters()->Rmax;

    // Invert z values
    //
    G4double* zValuesRefl = new G4double[nofZplanes];
    for (G4int i=0; i<nofZplanes; ++i)  { zValuesRefl[i] = - zValues[i]; }
    
    G4Polycone* newSolid
      = new G4Polycone(msol->GetName(),
                       msol->GetStartPhi(), 
                       msol->GetEndPhi() - msol->GetStartPhi(),
                       nofZplanes, zValuesRefl, rminValues, rmaxValues);

    delete [] zValuesRefl;       

    msol = newSolid;
    fmotherSolid = newSolid;
    fReflectedSolid = true;
    fDeleteSolid = true;
  }    
}

//---------------------------------------------------------------------
G4VParameterisationPolycone::~G4VParameterisationPolycone()
{
}

//---------------------------------------------------------------------
G4ParameterisationPolyconeRho::
G4ParameterisationPolyconeRho( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationPolycone( axis, nDiv, width, offset, msolid, divType )
{
  CheckParametersValidity();
  SetType( "DivisionPolyconeRho" );

  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( origparamMother->Rmax[0]
                         - origparamMother->Rmin[0], width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( origparamMother->Rmax[0]
                           - origparamMother->Rmin[0], nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= -1 )
  {
    G4cout << " G4ParameterisationPolyconeRho - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//---------------------------------------------------------------------
G4ParameterisationPolyconeRho::~G4ParameterisationPolyconeRho()
{
}

//---------------------------------------------------------------------
void G4ParameterisationPolyconeRho::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();

  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  if( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH )
  {
    std::ostringstream message;
    message << "In solid " << msol->GetName() << G4endl
            << "Division along R will be done with a width "
            << "different for each solid section." << G4endl
            << "WIDTH will not be used !";
    G4Exception("G4VParameterisationPolycone::CheckParametersValidity()",
                "GeomDiv1001", JustWarning, message);
  }
  if( foffset != 0. )
  {
    std::ostringstream message;
    message << "In solid " << msol->GetName() << G4endl
            << "Division along  R will be done with a width "
            << "different for each solid section." << G4endl
            << "OFFSET will not be used !";
    G4Exception("G4VParameterisationPolycone::CheckParametersValidity()",
                "GeomDiv1001", JustWarning, message);
  }
}

//------------------------------------------------------------------------
G4double G4ParameterisationPolyconeRho::GetMaxParameter() const
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
  G4PolyconeHistorical* original_pars = msol->GetOriginalParameters();
  return original_pars->Rmax[0] - original_pars->Rmin[0];
}


//---------------------------------------------------------------------
void
G4ParameterisationPolyconeRho::
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
    G4cout << " G4ParameterisationPolyconeRho " << G4endl
           << " foffset: " << foffset
           << " - fwidth: " << fwidth << G4endl;
  }
#endif

  ChangeRotMatrix( physVol );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconeRho "
           << G4endl
           << " Position: " << origin/mm
           << " - Width: " << fwidth/deg
           << " - Axis: " << faxis  << G4endl;
  }
#endif
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyconeRho::
ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyconeHistorical origparam( *origparamMother );
  G4int nZplanes = origparamMother->Num_z_planes;

  G4double width = 0.;
  for( G4int ii = 0; ii < nZplanes; ++ii )
  {
    width = CalculateWidth( origparamMother->Rmax[ii]
                          - origparamMother->Rmin[ii], fnDiv, foffset );
    origparam.Rmin[ii] = origparamMother->Rmin[ii]+foffset+width*copyNo;
    origparam.Rmax[ii] = origparamMother->Rmin[ii]+foffset+width*(copyNo+1);
  }

  pcone.SetOriginalParameters(&origparam);  // copy values & transfer pointers
  pcone.Reset();                            // reset to new solid parameters

#ifdef G4DIVDEBUG
  if( verbose >= -2 )
  {
    G4cout << "G4ParameterisationPolyconeRho::ComputeDimensions()" << G4endl
           << "-- Parametrised pcone copy-number: " << copyNo << G4endl;
    pcone.DumpInfo();
  }
#endif
}

//---------------------------------------------------------------------
G4ParameterisationPolyconePhi::
G4ParameterisationPolyconePhi( EAxis axis, G4int nDiv,
                               G4double width, G4double offset,
                               G4VSolid* msolid, DivisionType divType )
  :  G4VParameterisationPolycone( axis, nDiv, width, offset, msolid, divType )
{ 
  CheckParametersValidity();
  SetType( "DivisionPolyconePhi" );

  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
  G4double deltaPhi = msol->GetEndPhi() - msol->GetStartPhi();

  if( divType == DivWIDTH )
  {
    fnDiv = CalculateNDiv( deltaPhi, width, offset );
  }
  else if( divType == DivNDIV )
  {
    fwidth = CalculateWidth( deltaPhi, nDiv, offset );
  }

#ifdef G4DIVDEBUG
  if( verbose >= 1 )
  {
    G4cout << " G4ParameterisationPolyconePhi - # divisions " << fnDiv
           << " = " << nDiv << G4endl
           << " Offset " << foffset/deg << " = " << offset/deg << G4endl
           << " Width " << fwidth/deg << " = " << width/deg << G4endl;
  }
#endif
}

//---------------------------------------------------------------------
G4ParameterisationPolyconePhi::~G4ParameterisationPolyconePhi()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationPolyconePhi::GetMaxParameter() const
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);
  return msol->GetEndPhi() - msol->GetStartPhi();
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyconePhi::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  //----- translation 
  G4ThreeVector origin(0.,0.,0.); 
  //----- set translation 
  physVol->SetTranslation( origin );

  //----- calculate rotation matrix (so that all volumes point to the centre)
  G4double posi = foffset + copyNo*fwidth;

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyconePhi - position: " << posi/deg
           << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
#endif

  ChangeRotMatrix( physVol, -posi );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconePhi "
           << copyNo << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
#endif
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyconePhi::
ComputeDimensions( G4Polycone& pcone, const G4int,
                   const G4VPhysicalVolume* ) const
{
  G4Polycone* msol = (G4Polycone*)(fmotherSolid);

  G4PolyconeHistorical* origparamMother = msol->GetOriginalParameters();
  G4PolyconeHistorical origparam( *origparamMother );
  origparam.Start_angle = origparamMother->Start_angle;
  origparam.Opening_angle = fwidth;

  pcone.SetOriginalParameters(&origparam);  // copy values & transfer pointers
  pcone.Reset();                            // reset to new solid parameters

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << "G4ParameterisationPolyconePhi::ComputeDimensions():" << G4endl;
    pcone.DumpInfo();
  }
#endif
}

//---------------------------------------------------------------------
G4ParameterisationPolyconeZ::
G4ParameterisationPolyconeZ( EAxis axis, G4int nDiv,
                             G4double width, G4double offset,
                             G4VSolid* msolid, DivisionType divType)
  : G4VParameterisationPolycone( axis, nDiv, width, offset, msolid, divType ),
    fOrigParamMother(((G4Polycone*)fmotherSolid)->GetOriginalParameters())
{

  CheckParametersValidity();
  SetType( "DivisionPolyconeZ" );

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
    G4cout << " G4ParameterisationPolyconeZ - # divisions " << fnDiv << " = "
           << nDiv << G4endl
           << " Offset " << foffset << " = " << offset << G4endl
           << " Width " << fwidth << " = " << width << G4endl;
  }
#endif
}

//---------------------------------------------------------------------
G4ParameterisationPolyconeZ::~G4ParameterisationPolyconeZ()
{
}

//------------------------------------------------------------------------
G4double G4ParameterisationPolyconeZ::GetR(G4double z, 
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
G4double G4ParameterisationPolyconeZ::GetRmin(G4double z, G4int nseg) const
{
// Get Rmin in the given z position for the given polycone segment 

  return GetR(z, 
              fOrigParamMother->Z_values[nseg], 
              fOrigParamMother->Rmin[nseg],
              fOrigParamMother->Z_values[nseg+1], 
              fOrigParamMother->Rmin[nseg+1]);
}  
                                           
//------------------------------------------------------------------------
G4double G4ParameterisationPolyconeZ::GetRmax(G4double z, G4int nseg) const
{
// Get Rmax in the given z position for the given polycone segment 

  return GetR(z, 
              fOrigParamMother->Z_values[nseg], 
              fOrigParamMother->Rmax[nseg],
              fOrigParamMother->Z_values[nseg+1], 
              fOrigParamMother->Rmax[nseg+1]);
}  
                                           
//------------------------------------------------------------------------
G4double G4ParameterisationPolyconeZ::GetMaxParameter() const
{
  return std::abs (fOrigParamMother->Z_values[fOrigParamMother->Num_z_planes-1]
                  -fOrigParamMother->Z_values[0]);
}

//---------------------------------------------------------------------
void G4ParameterisationPolyconeZ::CheckParametersValidity()
{
  G4VDivisionParameterisation::CheckParametersValidity();
  
  // Division will be following the mother polycone segments
  //
  if( fDivisionType == DivNDIV )
  {
    if( fnDiv > fOrigParamMother->Num_z_planes-1 )
    { 
      std::ostringstream error;
      error  << "Configuration not supported." << G4endl
             << "Division along Z will be done by splitting in the defined"
             << G4endl
             << "Z planes, i.e, the number of division would be: "
             << fOrigParamMother->Num_z_planes-1
             << ", instead of: " << fnDiv << " !"; 
      G4Exception("G4ParameterisationPolyconeZ::CheckParametersValidity()",
                  "GeomDiv0001", FatalException, error);
    }
  }  
     
  // Division will be done within one polycone segment
  // with applying given width and offset
  //
  if( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH )
  {
    // Check if divided region does not span over more
    // than one z segment
  
    G4int isegstart = -1;  // number of the segment containing start position
    G4int isegend = -1;    // number of the segment containing end position

    if ( !fReflectedSolid )
    {
      // The start/end position of the divided region
      //
      G4double zstart 
        = fOrigParamMother->Z_values[0] + foffset;
      G4double zend 
        = fOrigParamMother->Z_values[0] + foffset + fnDiv* fwidth;
   
      G4int counter = 0;
      while ( isegend < 0 && counter < fOrigParamMother->Num_z_planes-1 )
      {
        // first segment
        if ( zstart >= fOrigParamMother->Z_values[counter]  &&
             zstart  < fOrigParamMother->Z_values[counter+1] )
        {
           isegstart = counter;
        }     
        // last segment
        if ( zend  > fOrigParamMother->Z_values[counter] &&
             zend <= fOrigParamMother->Z_values[counter+1] )
        {
          isegend = counter;
        }   
        ++counter;   
      }  // Loop checking, 06.08.2015, G.Cosmo
    }
    else
    {
      // The start/end position of the divided region
      //
      G4double zstart 
        = fOrigParamMother->Z_values[0] - foffset;
      G4double zend 
        = fOrigParamMother->Z_values[0] - ( foffset + fnDiv* fwidth);
   
      G4int counter = 0;
      while ( isegend < 0 && counter < fOrigParamMother->Num_z_planes-1 )
      {
        // first segment
        if ( zstart <= fOrigParamMother->Z_values[counter]  &&
             zstart  > fOrigParamMother->Z_values[counter+1] )
        {
           isegstart = counter;
        }     
        // last segment
        if ( zend  < fOrigParamMother->Z_values[counter] &&
             zend >= fOrigParamMother->Z_values[counter+1] )
        {
           isegend = counter;
        }   
        ++counter;   
      }  // Loop checking, 06.08.2015, G.Cosmo
    }
      
  
    if ( isegstart != isegend )
    {
      std::ostringstream message;
      message << "Condiguration not supported." << G4endl
              << "Division with user defined width." << G4endl
              << "Solid " << fmotherSolid->GetName() << G4endl
              << "Divided region is not between two z planes.";
      G4Exception("G4ParameterisationPolyconeZ::CheckParametersValidity()",
                  "GeomDiv0001", FatalException, message);
    }
  
    fNSegment = isegstart;
  }  
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyconeZ::
ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  if ( fDivisionType == DivNDIV )
  {
    // The position of the centre of copyNo-th mother polycone segment
    //
    G4double posi = ( fOrigParamMother->Z_values[copyNo]
                    + fOrigParamMother->Z_values[copyNo+1])/2;
    physVol->SetTranslation( G4ThreeVector(0, 0, posi) );
  }
  
  if ( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH )
  {
    // The position of the centre of copyNo-th division
    //
    G4double posi = fOrigParamMother->Z_values[0];
      
    if ( !fReflectedSolid )  
      posi += foffset + (2*copyNo + 1) * fwidth/2.;
    else
      posi -= foffset + (2*copyNo + 1) * fwidth/2.;
    
    physVol->SetTranslation( G4ThreeVector(0, 0, posi) );
  }   

  //----- calculate rotation matrix: unit

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << " G4ParameterisationPolyconeZ - position: " << posi << G4endl
           << " copyNo: " << copyNo << " - foffset: " << foffset/deg
           << " - fwidth: " << fwidth/deg << G4endl;
  }
#endif

  ChangeRotMatrix( physVol );

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << std::setprecision(8) << " G4ParameterisationPolyconeZ "
           << copyNo << G4endl
           << " Position: " << origin << " - Width: " << fwidth
           << " - Axis: " << faxis  << G4endl;
  }
#endif
}

//---------------------------------------------------------------------
void
G4ParameterisationPolyconeZ::
ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                   const G4VPhysicalVolume* ) const
{

  // Define division solid
  //
  G4PolyconeHistorical origparam;
  G4int nz = 2; 
  origparam.Num_z_planes = nz;
  origparam.Start_angle = fOrigParamMother->Start_angle;
  origparam.Opening_angle = fOrigParamMother->Opening_angle;

  // Define division solid z sections
  //
  origparam.Z_values = new G4double[nz];
  origparam.Rmin = new G4double[nz];
  origparam.Rmax = new G4double[nz];

  if ( fDivisionType == DivNDIV )
  {
    // The position of the centre of copyNo-th mother polycone segment
    G4double posi = (fOrigParamMother->Z_values[copyNo]
                   + fOrigParamMother->Z_values[copyNo+1])/2;

    origparam.Z_values[0] = fOrigParamMother->Z_values[copyNo] - posi;
    origparam.Z_values[1] = fOrigParamMother->Z_values[copyNo+1] - posi;
    origparam.Rmin[0] = fOrigParamMother->Rmin[copyNo];
    origparam.Rmin[1] = fOrigParamMother->Rmin[copyNo+1];
    origparam.Rmax[0] = fOrigParamMother->Rmax[copyNo];
    origparam.Rmax[1] = fOrigParamMother->Rmax[copyNo+1];
  }
  
  if ( fDivisionType == DivNDIVandWIDTH || fDivisionType == DivWIDTH )
  {
    if ( !fReflectedSolid )
    {
      origparam.Z_values[0] = - fwidth/2.;
      origparam.Z_values[1] = fwidth/2.;

      // The position of the centre of copyNo-th division
      //
      G4double posi = fOrigParamMother->Z_values[0]
                    + foffset + (2*copyNo + 1) * fwidth/2.;
    
      // The first and last z sides z values
      //
      G4double zstart = posi - fwidth/2.;
      G4double zend = posi + fwidth/2.;
      origparam.Rmin[0] = GetRmin(zstart, fNSegment); 
      origparam.Rmax[0] = GetRmax(zstart, fNSegment);  
      origparam.Rmin[1] = GetRmin(zend, fNSegment); 
      origparam.Rmax[1] = GetRmax(zend, fNSegment);  
    }
    else
    {
      origparam.Z_values[0] = fwidth/2.;
      origparam.Z_values[1] = - fwidth/2.;

      // The position of the centre of copyNo-th division
      //
      G4double posi = fOrigParamMother->Z_values[0]
                    - ( foffset + (2*copyNo + 1) * fwidth/2.);
    
      // The first and last z sides z values
      //
      G4double zstart = posi + fwidth/2.;
      G4double zend = posi - fwidth/2.;
      origparam.Rmin[0] = GetRmin(zstart, fNSegment); 
      origparam.Rmax[0] = GetRmax(zstart, fNSegment);  
      origparam.Rmin[1] = GetRmin(zend, fNSegment); 
      origparam.Rmax[1] = GetRmax(zend, fNSegment);  
    }

    // It can happen due to rounding errors
    //
    if ( origparam.Rmin[0]    < 0.0 ) origparam.Rmin[0] = 0.0;
    if ( origparam.Rmin[nz-1] < 0.0 ) origparam.Rmin[1] = 0.0;
  }  

  pcone.SetOriginalParameters(&origparam);  // copy values & transfer pointers
  pcone.Reset();                            // reset to new solid parameters

#ifdef G4DIVDEBUG
  if( verbose >= 2 )
  {
    G4cout << "G4ParameterisationPolyconeZ::ComputeDimensions()" << G4endl
           << "-- Parametrised pcone copy-number: " << copyNo << G4endl;
    pcone.DumpInfo();
  }
#endif
}
