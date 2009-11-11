/////////////////////////////////////////////////////////////////////////////////
//      Class:		G4AdjointCSMatrix.hh
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
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

        //////////////
        // Methods  // 
	//////////////
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
	std::vector< std::vector< size_t >* > theLogProbMatrixIndex; //index of equidistant LogProb
	std::vector< G4double > log0Vector;
	
	unsigned int nb_of_PrimEnergy;
	G4bool is_scat_proj_to_proj_case;
	G4double dlog;
	

};
#endif
