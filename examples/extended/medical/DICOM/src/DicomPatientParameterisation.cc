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
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Universit.AŽé Laval, QuŽébec (QC) Canada
//*******************************************************

#include "DicomPatientParameterisation.hh"
#include "DicomConfiguration.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

DicomPatientParameterisation::DicomPatientParameterisation(G4int, // NoVoxels, 
							   G4double maxDensity, 
							   G4double minDensity ,
							   G4Material* lunginhale,
							   G4Material* lungexhale,
							   G4Material* adipose,
							   G4Material* breast,
							   G4Material* phantom,
							   G4Material* muscle,
							   G4Material* liver,
							   G4Material* denseBone,
							   G4Material* trabecularBone) { 

    lungExhale    = lungexhale;
    lungInhale    = lunginhale;
    adiposeTissue = adipose;
    breastTissue  = breast;
    phantomTissue = phantom;
    muscleTissue  = muscle;
    liverTissue   = liver;
    denseBoneTissue = denseBone;
    trabecularBoneTissue = trabecularBone;
 
    density.clear();

    DicomConfiguration dicomConfiguration;

    // images must have the same dimension
    G4int totalNumberOfFile = dicomConfiguration.GetTotalNumberOfFile();
 
    G4double maxsl = -999. , minsl = 999.;
    for( G4int i = 0; i < totalNumberOfFile; i++) {
	G4double sliceLocation = dicomConfiguration.GetSliceLocation()[i];
	if(sliceLocation > maxsl) maxsl = sliceLocation;
	if(sliceLocation < minsl) minsl = sliceLocation;
    }
    middleLocationValue = (maxsl + minsl)*0.5;  

    //
    readColorChart();


    G4double alpha;

    attributeLungINhale = new G4VisAttributes;
    attributeLungINhale->SetColour(getChartColor(ColorChart::CRED, 0.217),
				   getChartColor(ColorChart::CGREEN, 0.217),
				   getChartColor(ColorChart::CBLUE, 0.217),
				   alpha=1.);
    attributeLungINhale->SetForceSolid(true);

    attributeLungEXhale = new G4VisAttributes;
    attributeLungEXhale->SetColour(getChartColor(ColorChart::CRED, 0.508),
				   getChartColor(ColorChart::CGREEN, 0.508),
				   getChartColor(ColorChart::CBLUE, 0.508),
				   alpha=1.);
    attributeLungEXhale->SetForceSolid(true);

    attributeAdipose = new G4VisAttributes;
    attributeAdipose->SetColour(getChartColor(ColorChart::CRED, 0.967),
				getChartColor(ColorChart::CGREEN, 0.967),
				getChartColor(ColorChart::CBLUE, 0.967),
				alpha=1.);
    attributeAdipose->SetForceSolid(true);

    attributeBreast = new G4VisAttributes;
    attributeBreast->SetColour(getChartColor(ColorChart::CRED, 0.99),
			       getChartColor(ColorChart::CGREEN, 0.99),
			       getChartColor(ColorChart::CBLUE, 0.99),
			       alpha=1.);
    attributeBreast->SetForceSolid(true);

    attributePhantom = new G4VisAttributes;
    attributePhantom->SetColour(getChartColor(ColorChart::CRED, 1.018),
				getChartColor(ColorChart::CGREEN, 1.018),
				getChartColor(ColorChart::CBLUE, 1.018),
				alpha=1.);
    attributePhantom->SetForceSolid(true);

    attributeMuscle = new G4VisAttributes;
    attributeMuscle->SetColour(getChartColor(ColorChart::CRED, 1.061),
			       getChartColor(ColorChart::CGREEN, 1.061),
			       getChartColor(ColorChart::CBLUE, 1.061),
			       alpha=1.);
    attributeMuscle->SetForceSolid(true);


    attributeLiver = new G4VisAttributes;
    attributeLiver->SetColour(getChartColor(ColorChart::CRED, 1.071),
			      getChartColor(ColorChart::CGREEN, 1.071),
			      getChartColor(ColorChart::CBLUE, 1.071),
			      alpha=1.);
    attributeLiver->SetForceSolid(true);

    attributeTrabecularBone = new G4VisAttributes;
    attributeTrabecularBone->SetColour(getChartColor(ColorChart::CRED, 1.159),
				       getChartColor(ColorChart::CGREEN, 1.159),
				       getChartColor(ColorChart::CBLUE, 1.159),
				       alpha=1.);
    attributeTrabecularBone->SetForceSolid(true);

    attributeDenseBone = new G4VisAttributes;
    attributeDenseBone->SetColour(getChartColor(ColorChart::CRED, 1.575),
				  getChartColor(ColorChart::CGREEN, 1.575),
				  getChartColor(ColorChart::CBLUE, 1.575),
				  alpha=1.);
    attributeDenseBone->SetForceSolid(true);

    attributeAir = new G4VisAttributes;
    attributeAir->SetColour(getChartColor(ColorChart::CRED, 1.),
			    getChartColor(ColorChart::CGREEN, 1.),
			    getChartColor(ColorChart::CBLUE, 1.),
			    alpha=1.);
    attributeAir->SetForceSolid(false);
  

    rows           = dicomConfiguration.GetTotalRows();
    columns        = dicomConfiguration.GetTotalColumns();
    compression    = dicomConfiguration.GetCompressionValue();
    max            = dicomConfiguration.GetTotalNumberOfFile();
    pixelSpacingX  = dicomConfiguration.GetXPixelSpacing();
    pixelSpacingY  = dicomConfiguration.GetYPixelSpacing();
    sliceThickness = dicomConfiguration.GetSliceThickness();

    GetDensity( maxDensity , minDensity );

    dicomConfiguration.ClearDensityData();
}

DicomPatientParameterisation::~DicomPatientParameterisation() {

    // visualisation attributes ...
    delete attributeAdipose;
    delete attributeLungEXhale;
    delete attributeBreast;
    delete attributePhantom;
    delete attributeMuscle;
    delete attributeLiver;
    delete attributeTrabecularBone;
    delete attributeLungINhale;
    delete attributeDenseBone;
    delete attributeAir;
  
    // materials ...
    delete trabecularBoneTissue;
    delete denseBoneTissue;
    delete liverTissue;
    delete muscleTissue;

    delete breastTissue; 
    delete adiposeTissue;
    delete lungInhale;
    delete lungExhale;
    
}

void DicomPatientParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
    G4double originZ = patientPlacementZ[copyNo]*mm-middleLocationValue*mm-sliceThickness/2;
    G4ThreeVector origin( patientPlacementX[copyNo]*mm, 
			  patientPlacementY[copyNo]*mm, 
			  originZ*mm );

    physVol->SetTranslation(origin);
}

void DicomPatientParameterisation::ComputeDimensions(G4Box& voxels, const G4int, const G4VPhysicalVolume*) const
{
    voxels.SetXHalfLength((pixelSpacingX * compression/2.0) * mm);
    voxels.SetYHalfLength((pixelSpacingY * compression/2.0) * mm);
    voxels.SetZHalfLength((sliceThickness / 2.0) * mm);
}

G4Material*
DicomPatientParameterisation::ComputeMaterial(const G4int copyNo,
                                              G4VPhysicalVolume* physVol,
                                              const G4VTouchable*)
{
    if( density[copyNo] >= 0.207 && density[copyNo] <= 0.227 ) {
	physVol->SetName("PhysicalLungINhale");
	physVol->GetLogicalVolume()->SetVisAttributes( attributeLungINhale );
	return lungInhale;

    } else if( density[copyNo] >= 0.481 && density[copyNo] <= 0.534 ) {
	physVol->SetName("PhysicalLungEXhale");
	physVol->GetLogicalVolume()->SetVisAttributes( attributeLungEXhale );
	return lungExhale;

    } else if( density[copyNo] >= 0.919 && density[copyNo] <= 0.979 ) {
	physVol->SetName("PhysicalAdipose");
	physVol->GetLogicalVolume()->SetVisAttributes( attributeAdipose );
	return adiposeTissue;

    } else if( density[copyNo] > 0.979 && density[copyNo] <= 1.004 ) {
	physVol->SetName("PhysicalBreast");
	physVol->GetLogicalVolume()->SetVisAttributes( attributeBreast );
	return breastTissue;

    } else if( density[copyNo] > 1.004 && density[copyNo] <= 1.043 ) {
	physVol->SetName("PhysicalPhantom");
	physVol->GetLogicalVolume()->SetVisAttributes( attributePhantom );
	return phantomTissue;

    } else if( density[copyNo] > 1.043 && density[copyNo] <= 1.109 ) {
	physVol->SetName("PhysicalMuscle");
	physVol->GetLogicalVolume()->SetVisAttributes( attributeMuscle );
	return muscleTissue;

    } else if( density[copyNo] > 1.109 && density[copyNo] <= 1.113 ) {
	physVol->SetName("PhysicalLiver");
	physVol->GetLogicalVolume()->SetVisAttributes( attributeLiver );
	return liverTissue;

    } else if( density[copyNo] > 1.113 && density[copyNo] <= 1.217 ) {
	physVol->SetName("PhysicalTrabecularBone");
	physVol->GetLogicalVolume()->SetVisAttributes( attributeTrabecularBone );
	return trabecularBoneTissue;

    } else if( density[copyNo] > 1.496 && density[copyNo] <= 1.654 ) {
	physVol->SetName("PhysicalDenseBone");
	physVol->GetLogicalVolume()->SetVisAttributes( attributeDenseBone );
	return denseBoneTissue;

    }

    return physVol->GetLogicalVolume()->GetMaterial();
}

void DicomPatientParameterisation::GetDensity(G4double maxdensity, G4double mindensity) {

    DicomConfiguration dicomConfiguration;
  
    G4int copyCounter = 0;

    G4int totalNumberOfFile = dicomConfiguration.GetTotalNumberOfFile();

    G4int lenRows       = rows/compression;
    G4int lenColumns    = columns/compression;
    G4double xDimension = (lenColumns*pixelSpacingX)/2;

    G4int i = 0;
    for( G4int z = 0; z < totalNumberOfFile; z++ ) {

	G4double slicePosition = dicomConfiguration.GetSliceLocation()[z];      
	
	for( G4int j = 1; j <= lenRows; j++ ) {
	    for( G4int w = 1; w <= lenColumns; w++ ) {
	      
		G4double tissueDensity = dicomConfiguration.GetDensityValue(i++);

		if( tissueDensity != -1 ) {
		    if( tissueDensity >= mindensity && tissueDensity <= maxdensity ) {

			density.push_back( tissueDensity );
			copyCounter++;
			G4double yPixel = (pixelSpacingY/2 + (w-1)*pixelSpacingY);
			G4double yDimension = ((lenRows*pixelSpacingX)/2)-(pixelSpacingY/2+(j-1)*pixelSpacingY);
                      
			patientPlacementX.push_back( ( compression*(xDimension- yPixel ) ) *mm );
			patientPlacementY.push_back( ( compression* yDimension  ) *mm );
			patientPlacementZ.push_back( ( slicePosition + sliceThickness/2 ) *mm );
		    }
		}
	    }            
	}
    }

}

void DicomPatientParameterisation::readColorChart() {
    std::ifstream cm("Colormap.dat");
    if(!cm) {
	G4cerr << "Colormap.dat couldn't be opened!!" << G4endl;
	numColorChart = 2;
	ColorChart cc;
	cc.density = 0.;
	cc.color[ColorChart::CRED] = 0.;
	cc.color[ColorChart::CGREEN] = 0.;
	cc.color[ColorChart::CBLUE] = 0.;
	cc.alpha = 1.;
	colorChart.push_back(cc);
	cc.density = 4.;
	cc.color[ColorChart::CRED] = 1.;
	cc.color[ColorChart::CGREEN] = 1.;
	cc.color[ColorChart::CBLUE] = 1.;
	cc.alpha = 1.;
	colorChart.push_back(cc);
	return;
    }

    cm >> numColorChart;
    ColorChart cc;
    for(int i = 0; i < numColorChart; i++) {
	cm >> cc.density
	   >> cc.color[ColorChart::CRED]
	   >> cc.color[ColorChart::CGREEN]
	   >> cc.color[ColorChart::CBLUE]
	   >> cc.alpha; 
	colorChart.push_back(cc);
    }

}

G4double DicomPatientParameterisation::getChartColor(G4int CC, G4double density) {
    G4double color = 0.;

    for(int i = 0; i < numColorChart; i++) {
	if(density <= colorChart[i].density) {

	    G4double w = (density - colorChart[i-1].density)
		/(colorChart[i].density - colorChart[i-1].density);
	    color = w*colorChart[i-1].color[CC] + (1-w)*colorChart[i].color[CC];
	    return color;
	}
    }
    return color;
}
