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
/// \file SDDData.cc
/// \brief Implementation of the SDDData class 

#include "SDDData.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

SDDData::SDDData(std::string p_name):
filename_(p_name)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDDData::SDDHeader SDDData::ReadHeader()
{

	std::ifstream file(filename_.c_str());
	
	if(!file.is_open())
	{
		std::cout << "No file: " << filename_ << std::endl;
		return header_;
	}
	
	//1-SDD version
	ReadString(file,header_.sdd_version);
	//2-Software
	ReadString(file,header_.software);
	//3-Author
	ReadString(file,header_.author);
	//4-Simulation details
	ReadString(file,header_.sim_details);

	//5- Source
	ReadString(file,header_.src_details);
	//6-Source type
	ReadInt(file,header_.src_type);
	//7-Incident particles
	ReadInts(file,header_.src_pdg);
	//8-Mean Particle energy
	ReadDoubles(file,header_.src_energy);
	//9-Energy distribution
	ReadString(file,header_.energy_dist);
	//10-Particle fraction
	ReadDoubles(file,header_.part_fraction);
	//11-Dose or fluence
	ReadDoubles(file,header_.dose);
	//12-Dose rate
	ReadDouble(file,header_.dose_rate);

	//13-Irradiation target
	ReadString(file,header_.target);
	//14-Volumes
	ReadDoubles(file,header_.volumes);
	//15-Chromosome sizes
	ReadDoubles(file,header_.chromo_size);
	//16-DNA density
	ReadDouble(file,header_.dna_density);

	//17-Cell cycle phase
	ReadDoubles(file,header_.cell_cycle);

	//18-DNA strcuture
	ReadInts(file,header_.dna_struct);

	//19- in vitro/in vivo
	ReadInt(file,header_.vitro_vivo);
	
	//20-Proliferation status
	ReadString(file,header_.proliferation);
	
	//21-Microenvironment
	ReadDoubles(file,header_.microenv);

	//22-Damage definition
	ReadDoubles(file,header_.damage_def);
	
	//23-Time
	ReadDouble(file,header_.time);
	
	//24-Damage and primary count
	ReadInts(file,header_.damage_prim_count);
	
	//25-Data entries
	ReadBools(file,header_.entries);

	//26-Additional information
	ReadString(file,header_.info);
	
	file.close();

	return header_;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDData::ParseData()
{
	ReadHeader();
	
	data_.clear();

	std::ifstream file(filename_.c_str());
	
	std::string line;

	// Pass the header
	while(std::getline(file,line))
	{
		if(line.find("EndOfHeader")!=std::string::npos)
			break;
	}

	// Start to read the data

	while(std::getline(file,line))
	{
		if(!line.empty())
			ParseLineData(line);
	}
	
	file.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDData::ParseLineData(std::string& line)
{

	SDDDamage dmg;

	if(header_.entries[0])
		ExtractInts(line,2,dmg.classification);
	if(header_.entries[1])
		ExtractDoubles(line,3,dmg.coordinates);
	if(header_.entries[2])
		ExtractInts(line,4,dmg.chromo_ID);
	if(header_.entries[3])
		ExtractDoubles(line,1,dmg.chromo_position);
	if(header_.entries[4])
		ExtractInts(line,3,dmg.cause);
	if(header_.entries[5])
		ExtractInts(line,3,dmg.types);
	if(header_.entries[6])
		ExtractInts(line,3,dmg.break_spec);
	if(header_.entries[7])
		ExtractInts(line,3,dmg.dna_seq);
	if(header_.entries[8])
		ExtractDoubles(line,1,dmg.lesion_time);
	if(header_.entries[9])
		ExtractInts(line,1,dmg.particles);
	if(header_.entries[10])
		ExtractDoubles(line,1,dmg.energy);
	if(header_.entries[11])
		ExtractDoubles(line,3,dmg.translation);
	if(header_.entries[12])
		ExtractDoubles(line,1,dmg.direction);
	if(header_.entries[13])
		ExtractDoubles(line,1,dmg.particle_time);

	data_.push_back(dmg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > SDDData::GetAllDamage()
{
	if(data_.size()<=0)
	{
		ParseData();
	}
	
	std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > fmDamage;
	
	for(auto it=data_.begin();it!=data_.end();it++)
	{

		Damage::DamageType pType = Damage::DamageType::fOther;
		unsigned int pChromo = -1;
		unsigned int pEvt = -1;
		unsigned int pStrand = -1;
		unsigned long int pCopyNb = -1;
		Position pPos(0,0,0);
		Damage::DamageCause pCause = Damage::DamageCause::fUnknown;
		Damage::DamageChromatin pChromatin = Damage::DamageChromatin::fUnspecified;

		if(header_.entries[0])
		{
			pEvt = it->classification[1];
		}
		if(header_.entries[1])
		{
			pPos.setX(it->coordinates[0]);
			pPos.setY(it->coordinates[1]);
			pPos.setZ(it->coordinates[2]);
		}
		if(header_.entries[2])
		{
			switch(it->chromo_ID[0])
			{
				case 0:
				pChromatin = Damage::DamageChromatin::fUnspecified;
				break;
				case 1:
				pChromatin = Damage::DamageChromatin::fHeterochromatin;
				break;
				case 2:
				pChromatin = Damage::DamageChromatin::fEuchromatin;
				break;
				case 3:
				pChromatin = Damage::DamageChromatin::fFreeDNA;
				break;
				case 4:
				pChromatin = Damage::DamageChromatin::fOtherDNA;
				break;
				default:
				pChromatin = Damage::DamageChromatin::fUnspecified;
			}
			pChromo = it->chromo_ID[1];
			pStrand = it->chromo_ID[3];
		}
		if(header_.entries[3])
		{
			pCopyNb = (unsigned long int)(it->chromo_position[0]);
		}
		if(header_.entries[4])
		{
			switch(it->cause[0])
			{
				case 0:
				pCause = Damage::DamageCause::fDirect;
				break;
				case 1:
				pCause = Damage::DamageCause::fIndirect;
				break;
				default:
				pCause = Damage::DamageCause::fUnknown;

			}
		}
		if(header_.entries[5])
		{
			if(it->types[0]>0)
				pType = Damage::DamageType::fBase;
			if(it->types[1]>0)
				pType = Damage::DamageType::fBackbone;
		}
		Damage aDamage(pType,pChromo,pEvt,pStrand,pCopyNb,pPos,pCause,pChromatin);
		auto chroPos = fmDamage.find(pChromo);
		if (chroPos == fmDamage.end()) {
			std::vector<Damage> dmv{aDamage};
			std::map<unsigned int,std::vector<Damage> > evtDamages{{pEvt,dmv}};
			fmDamage.insert({pChromo,evtDamages});
		} else {
			auto evtPos = chroPos->second.find(pEvt);
			if (evtPos == chroPos->second.end()) {
				std::vector<Damage> dmv{aDamage};
				chroPos->second.insert({pEvt,dmv});
			} else {
				chroPos->second[pEvt].push_back(aDamage);
			}
		}
	}
	
	return fmDamage;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDData::ExtractInts(std::string& strLine,int numInt,std::vector<int>& field)
{
	std::string delimiter = ";";
	std::string token;

	size_t pos = strLine.find(delimiter);

	if(pos!=std::string::npos)
	{
		token = strLine.substr(0,pos);
		strLine.erase(0,pos+delimiter.length());
	}
	else
	{
		token = strLine;
		strLine = "";
	}

	std::string value;
	delimiter = ",";

	for(int i=1;i<numInt;i++)
	{
		pos = token.find(delimiter);
		value = token.substr(0,pos);
		field.push_back(std::atoi(value.c_str()));
		token.erase(0,pos+delimiter.length());
	}

	field.push_back(std::atoi(token.c_str()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDData::ExtractDoubles(std::string& strLine,int numInt,std::vector<double>& field)
{
	std::string delimiter = ";";
	std::string token;

	size_t pos = strLine.find(delimiter);

	if(pos!=std::string::npos)
	{
		token = strLine.substr(0,pos);
		strLine.erase(0,pos+delimiter.length());
	}
	else
	{
		token = strLine;
		strLine = "";
	}

	std::string value;
	delimiter = ",";

	for(int i=1;i<numInt;i++)
	{
		pos = token.find(delimiter);
		value = token.substr(0,pos);
		field.push_back(std::atoi(value.c_str()));
		token.erase(0,pos+delimiter.length());
	}

	field.push_back(std::stod(token.c_str()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDData::ReadString(std::ifstream& file,std::string& field)
{
	std::string line;
	std::string token;

	std::getline(file,line);
	std::istringstream ss(line);

	std::getline(ss,token,',');
	std::getline(ss,token,',');

	field=token;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDData::ReadInt(std::ifstream& file,int& field)
{
	std::string line;
	std::string token;

	std::getline(file,line);
	line = line.substr(0, line.size()-1);

	std::istringstream ss(line);

	std::getline(ss,token,',');
	std::getline(ss,token,',');

	field=std::atoi(token.c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDData::ReadInts(std::ifstream& file,std::vector<int>& field)
{
	std::string line;
	std::string token;

	std::getline(file,line);
	line = line.substr(0, line.size()-1);

	std::istringstream ss(line);

	std::getline(ss,token,',');

	while(std::getline(ss,token,','))
		field.push_back(std::atoi(token.c_str()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDData::ReadDouble(std::ifstream& file,double& field)
{
	std::string line;
	std::string token;

	std::getline(file,line);
	line = line.substr(0, line.size()-1);

	std::istringstream ss(line);

	std::getline(ss,token,',');
	std::getline(ss,token,',');

	field=std::stod(token.c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDData::ReadDoubles(std::ifstream& file,std::vector<double>& field)
{
	std::string line;
	std::string token;

	std::getline(file,line);
	line = line.substr(0, line.size()-1);

	std::istringstream ss(line);

	std::getline(ss,token,',');

	while(std::getline(ss,token,','))
		field.push_back(std::stod(token.c_str()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDData::ReadBools(std::ifstream& file,std::vector<bool>& field)
{
	std::string line;
	std::string token;

	std::getline(file,line);
	line = line.substr(0, line.size()-1);

	std::istringstream ss(line);

	std::getline(ss,token,',');

	while(std::getline(ss,token,','))
		field.push_back(std::stoi(token.c_str()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double SDDData::GetDose()
{
	return header_.dose[1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::map<int,unsigned long long int> SDDData::GetChromosomeBpSizesMap(double &sum)
{
	sum=0;
	std::map<int,unsigned long long int> chromap;
	for (int i=1;i<header_.chromo_size.size(); i++) { // index 0 of header_.chromo_size is number of chromosomes
		double nbp = header_.chromo_size[i]; // in Mbp
		nbp *= 1e+6; // to bp
		sum += nbp;
		chromap.insert({(i-1),nbp}); 
	}
	return chromap;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
