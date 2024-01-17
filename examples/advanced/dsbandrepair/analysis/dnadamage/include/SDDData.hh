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
//
/// \file SDDData.hh
/// \brief Definition of the SDDData class

#ifndef SDDDATA_HH
#define SDDDATA_HH

#include <string>
#include <vector>
#include <sstream>
#include <map>
#include "Damage.hh"

class SDDData
{
public:

/** SDD header structure
 *
 *  Define the content of SSD file header
 */
	struct SDDHeader 
	{  
		std::string        sdd_version{""};

		std::string         software{""};

		std::string         author{""};

		std::string         sim_details{""};

		std::string         src_details{""};
		int                 src_type{0};
		std::vector<int>    src_pdg;
		std::vector<double> src_energy;
		std::string         energy_dist{""};
		std::vector<double> part_fraction;
		std::vector<double> dose;
		double              dose_rate{0};

		std::string         target{""};
		std::vector<double> volumes;
		std::vector<double> chromo_size;
		double              dna_density{0};

		std::vector<double> cell_cycle;
		std::vector<int>    dna_struct;
		int                 vitro_vivo{0};
		std::string         proliferation{""};
		std::vector<double> microenv;

		std::vector<double> damage_def;
		double              time{0};
		std::vector<int>    damage_prim_count;

		std::vector<bool>   entries;

		std::string         info{""};
	}; 

/** SDD damage structure
 *
 *  Define the content of a damage as stored in SDD file
 */
	struct SDDDamage
	{
		std::vector<int>		classification;
		std::vector<double>		coordinates;
		std::vector<int>		chromo_ID;
		std::vector<double>		chromo_position;
		std::vector<int>		cause;
		std::vector<int>		types;

		std::vector<int>		break_spec;
		std::vector<int>		dna_seq;
		std::vector<double>		lesion_time;

		std::vector<int>		particles;
		std::vector<double> 	energy;
		std::vector<double>		translation;
		std::vector<double>		direction;
		std::vector<double>		particle_time;
	}; 

  	/** Constructor */
 	SDDData(std::string /*p_name*/);
  	/** Destructor */
	~SDDData() = default;

	/** Read header
	*   Reads and returns the header of a SDD file
	*/
	SDDHeader ReadHeader();

	/** Chromosome sizes
	* Returns the list of sizes of each chromosome in the cell geometry
	*/
	std::map<int,unsigned long long int> GetChromosomeBpSizesMap(double &sum);

	/** Dose
	* Returns the absorbed dose
	*/
	double GetDose();

	/** Parse data
	*   Parse data of SDD files and stores the data
	*   as SDDDamage in data_
	*/
	void ParseData();
	/** Get SDD damage
	* 	Returns all the SDD damage (i.e. data_)
	*/
	std::vector<SDDDamage>& GetSDDDamage(){return data_;};

	/** Get Damage
	*   Returns all the damage that have been converted into lighter object
	*   see Damage class
	*/
	std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > GetAllDamage();

private:

	std::string filename_;
	SDDHeader header_;
	std::vector<SDDDamage> data_;

	void ParseLineData(std::string&);

	void ReadString(std::ifstream&,std::string&);
	void ReadInt(std::ifstream&,int&);
	void ReadInts(std::ifstream&,std::vector<int>&);
	void ReadDouble(std::ifstream&,double&);
	void ReadDoubles(std::ifstream&,std::vector<double>&);
	void ReadBools(std::ifstream&,std::vector<bool>&);

	void ExtractInts(std::string&,int,std::vector<int>&);
	void ExtractDoubles(std::string&,int,std::vector<double>&);

};

#endif // SDDDATA_HH
