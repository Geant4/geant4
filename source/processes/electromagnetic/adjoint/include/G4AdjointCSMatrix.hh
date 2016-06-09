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
/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointCSMatrix.hh
//	Author:       	L. Desorgher
//	Date:		1st April 2007
// 	Organisation: 	SpaceIT GmbH
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	1st April 2007 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		An adjoint CS matrix is used by the model of a reverse process to sample an adjoint secondary (being equivalent to a forward primary). 
//		It represents the integration over the energy of the adjoint secondary (therefore the forward primary) of the differential cross section 
//		of the equiavlent forward  discrete process (Ionisation, Brem, PE effect, Compton,..) . Each reverse model has its own cross section matrix for a given cut, 
//		material couple. It is therefore recompute after a modification  of the cuts by the user. 
//		
//		
//

#ifndef G4AdjointCSMatrix_h
#define G4AdjointCSMatrix_h 1

#include"globals.hh"
#include<vector>
#include"G4ParticleDefinition.hh"

////////////////////////////////////////////////////////////////////////////////
//
class G4AdjointCSMatrix
{
        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////
public:
	G4AdjointCSMatrix(G4bool aBool);
	~G4AdjointCSMatrix();

        ////////////
        // Methods
	////////////
	void Clear();
	void AddData(G4double aPrimEnergy,G4double aCS, std::vector< G4double>* aLogSecondEnergyVector,
	 					        std::vector< G4double>* aLogProbVector,size_t n_pro_decade=0);	
	
	bool GetData(unsigned int i, G4double& aPrimEnergy,G4double& aCS,G4double& log0, std::vector< G4double>*& aLogSecondEnergyVector,
	 							      std::vector< G4double>*& aLogProbVector,
								      std::vector< size_t>*& aLogProbVectorIndex);
	
	inline std::vector< G4double >* GetLogPrimEnergyVector(){return &theLogPrimEnergyVector;}
	inline std::vector< G4double >* GetLogCrossSectionvector(){return &theLogCrossSectionVector;}
	inline G4double GetDlog(){return dlog;} 	
	inline G4bool IsScatProjToProjCase(){return is_scat_proj_to_proj_case;} 
	void Write(G4String file_name);
	void Read(G4String file_name);		

private:
        
	// we did first try to use G4PhysicsOrderedVector but they are not general enough for our purpose
	
	std::vector< G4double > theLogPrimEnergyVector; 
        std::vector< G4double > theLogCrossSectionVector; //Adjoint Cross sections in function of primary energy
        std::vector< std::vector< G4double >* > theLogSecondEnergyMatrix;
	std::vector< std::vector< G4double >* > theLogProbMatrix; //Each column represents the integrated probability of getting a secondary 
								      // in function of their energy 
	std::vector< std::vector< size_t >* > theLogProbMatrixIndex; //index of euqidistant LogProb
	std::vector< G4double > log0Vector;
	
	unsigned int nb_of_PrimEnergy;
	G4bool is_scat_proj_to_proj_case;
	G4double dlog;
	

};
#endif
