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
// $Id: G4OpticalSurface.hh 104459 2017-05-31 15:54:52Z gcosmo $
//
// 
////////////////////////////////////////////////////////////////////////
// G4OpticalSurface Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpticalSurface.hh
// Description: A optical surface class for use in G4OpBoundaryProcess
// Version:     2.0
// Created:     1997-06-26
// Author:      Peter Gumplinger
// Updated:     1999-10-29 add method and class descriptors
//              2017-02-24 Mariele Stockhoff add DAVIS model 
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4OpticalSurface_h
#define G4OpticalSurface_h 1

/////////////
// Includes
/////////////

#include "G4Types.hh"
#include "G4Physics2DVector.hh"
#include "G4SurfaceProperty.hh"

// Class Description:
// A optical surface class for use in the G4OpBoundaryProcess class.
// Contains the enumerations: G4OpticalSurfaceFinish, G4OpticalSurfaceType,
// and G4OpticalSurfaceModel.
// Class Description - End:

enum G4OpticalSurfaceFinish
{
   polished,                    // smooth perfectly polished surface
   polishedfrontpainted,        // smooth top-layer (front) paint
   polishedbackpainted,         // same is 'polished' but with a back-paint

   ground,                      // rough surface
   groundfrontpainted,          // rough top-layer (front) paint
   groundbackpainted,           // same as 'ground' but with a back-paint

   polishedlumirrorair,         // mechanically polished surface, with lumirror
   polishedlumirrorglue,        // mechanically polished surface, with lumirror & meltmount
   polishedair,                 // mechanically polished surface
   polishedteflonair,           // mechanically polished surface, with teflon
   polishedtioair,              // mechanically polished surface, with tio paint
   polishedtyvekair,            // mechanically polished surface, with tyvek
   polishedvm2000air,           // mechanically polished surface, with esr film
   polishedvm2000glue,          // mechanically polished surface, with esr film & meltmount

   etchedlumirrorair,           // chemically etched surface, with lumirror
   etchedlumirrorglue,          // chemically etched surface, with lumirror & meltmount
   etchedair,                   // chemically etched surface
   etchedteflonair,             // chemically etched surface, with teflon
   etchedtioair,                // chemically etched surface, with tio paint
   etchedtyvekair,              // chemically etched surface, with tyvek
   etchedvm2000air,             // chemically etched surface, with esr film
   etchedvm2000glue,            // chemically etched surface, with esr film & meltmount

   groundlumirrorair,           // rough-cut surface, with lumirror
   groundlumirrorglue,          // rough-cut surface, with lumirror & meltmount
   groundair,                   // rough-cut surface
   groundteflonair,             // rough-cut surface, with teflon
   groundtioair,                // rough-cut surface, with tio paint
   groundtyvekair,              // rough-cut surface, with tyvek
   groundvm2000air,             // rough-cut surface, with esr film
   groundvm2000glue,            // rough-cut surface, with esr film & meltmount
   
   // for DAVIS model
   Rough_LUT,                   //rough surface 
   RoughTeflon_LUT,             //rough surface wrapped in Teflon tape
   RoughESR_LUT,                //rough surface wrapped with ESR
   RoughESRGrease_LUT,          //rough surface wrapped with ESR
                                //and coupled with opical grease
   Polished_LUT,                //polished surface 
   PolishedTeflon_LUT,          //polished surface wrapped in Teflon tape
   PolishedESR_LUT,             //polished surface wrapped with ESR
   PolishedESRGrease_LUT,       //polished surface wrapped with ESR
                                //and coupled with opical grease
   Detector_LUT                 //polished surface with optical grease 
};

enum G4OpticalSurfaceModel
{
   glisur,                      // original GEANT3 model
   unified,                     // UNIFIED model
   LUT,                         // Look-Up-Table model
   DAVIS,                       // DAVIS model
   dichroic                     // dichroic filter
};

class G4MaterialPropertiesTable;

/////////////////////
// Class Definition
/////////////////////

class G4OpticalSurface : public G4SurfaceProperty
{

public: // Without description
  
        //////////////
        // Operators
        //////////////
  
	G4OpticalSurface(const G4OpticalSurface &right);
	G4OpticalSurface & operator=(const G4OpticalSurface &right);
  
	G4int operator==(const G4OpticalSurface &right) const;
	G4int operator!=(const G4OpticalSurface &right) const;

public: // With description

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

	G4OpticalSurface(const G4String& name,
                         G4OpticalSurfaceModel model = glisur,
                         G4OpticalSurfaceFinish finish = polished,
                         G4SurfaceType type = dielectric_dielectric,
                         G4double value = 1.0);
        // Constructor of an optical surface object.

public: // Without description

	virtual ~G4OpticalSurface();

	////////////
	// Methods
        ////////////

	// public methods

public: // With description

        void SetType(const G4SurfaceType& type);

        inline G4OpticalSurfaceFinish GetFinish() const { return theFinish; }
        // Returns the optical surface finish.
        void SetFinish(const G4OpticalSurfaceFinish );
        // Sets the optical surface finish.

        inline G4OpticalSurfaceModel GetModel() const { return theModel; }
        // Returns the optical surface model used.
        inline void SetModel(const G4OpticalSurfaceModel model)
                                                      { theModel = model; }
        // Sets the optical surface model to be followed.

	inline G4double GetSigmaAlpha() const { return sigma_alpha; }
        // Returns an unified model surface parameter.
	inline void SetSigmaAlpha(const G4double s_a) { sigma_alpha = s_a; }
        // Sets an unified model surface parameter.

	G4double GetPolish() const { return polish; }
        // Returns the optical surface polish type.
	inline void SetPolish(const G4double plsh) { polish=plsh; }
        // Sets the optical surface polish type.

	inline G4MaterialPropertiesTable* GetMaterialPropertiesTable() const
                                       { return theMaterialPropertiesTable; }
        // Retrieves the pointer of the G4MaterialPropertiesTable 
        // attached to optical surface.

	inline void SetMaterialPropertiesTable(G4MaterialPropertiesTable *anMPT)
                                       { theMaterialPropertiesTable = anMPT; }
        // Attaches a G4MaterialPropertiesTable to the optical surface.

	void DumpInfo() const;
        // Prints information about the optical surface.

        void ReadLUTFile(void);
        // Method to read the Look-Up-Table into array AngularDistribution

        inline G4double GetAngularDistributionValue(G4int, G4int, G4int);
        
        // for DAVIS model

        inline G4double GetAngularDistributionValueLUT(G4int);
        // Returns the AngularDistributionValue

        void ReadLUTDAVISFile(void);
        // Method to read the Davis Look-Up-Table into array AngularDistribution
 
        void ReadReflectivityLUTFile(void);
        // Method to read the Look-Up-Table for reflectivity
        
        inline G4double GetReflectivityLUTValue(G4int);
        // Returns the reflectivity value from the Davis Look-Up-Table
        
        inline G4int GetInmax(void) const { return indexmax; }
        // Returns the number of lines in the Davis Look-Up-Table

        inline G4int GetLUTbins(void) const { return LUTbins; }
        // Returns the number of probability values per incidentangle

        inline G4int GetRefMax(void) const { return RefMax; }
        // Returns the number of reflectivity values per angle
    
        inline G4int GetThetaIndexMax(void) const { return thetaIndexMax; }
        inline G4int GetPhiIndexMax(void) const { return phiIndexMax; }

        void ReadDichroicFile(void);
        // Method to read the dichroic surface data file into Dichroic

        inline G4Physics2DVector* GetDichroicVector();

private:

// ------------------
// Basic data members ( To define an optical surface)
// ------------------

        G4OpticalSurfaceModel theModel;		// Surface model
        G4OpticalSurfaceFinish theFinish;	// Surface finish

	G4double sigma_alpha;		// The sigma of micro-facet polar angle
	G4double polish;		// Polish parameter in glisur model

	G4MaterialPropertiesTable* theMaterialPropertiesTable;

        static const G4int incidentIndexMax = 91;
        static const G4int thetaIndexMax = 45;
        static const G4int phiIndexMax = 37;

        G4float* AngularDistribution;
        G4Physics2DVector* DichroicVector;

        // for DAVIS model
        static const G4int indexmax = 7280001; // 3640001; 
        static const G4int RefMax = 90; 
        static const G4int LUTbins =20000;
        G4float* AngularDistributionLUT;
        G4float* Reflectivity;
};

////////////////////
// Inline methods
////////////////////

inline
 G4double G4OpticalSurface::GetAngularDistributionValue(G4int angleIncident,
                                                        G4int thetaIndex,
                                                        G4int phiIndex)
{
  return AngularDistribution[angleIncident+
                             thetaIndex*incidentIndexMax+
                             phiIndex*thetaIndexMax*incidentIndexMax];
}

inline
 G4double G4OpticalSurface::GetAngularDistributionValueLUT(G4int i)
{
  return AngularDistributionLUT[i];
}

inline
 G4double G4OpticalSurface::GetReflectivityLUTValue(G4int i)
{
  return Reflectivity[i];
}

inline
 G4Physics2DVector* G4OpticalSurface::GetDichroicVector()
{
  return DichroicVector;
}

#endif /* G4OpticalSurface_h */
