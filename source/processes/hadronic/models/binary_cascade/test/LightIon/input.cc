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
#include <iostream>
#include <vector>
#include <utility>
#include <strstream>
#include <fstream>

int main()
{

   std::vector<std::pair<int,int> > materials;
   std::vector<std::pair<int,int> > projectiles;
   std::vector<double> energies;
   
   std::pair<int,int> H(1,1);
   std::pair<int,int> d(2,1);
   std::pair<int,int> O(16,8);
   std::pair<int,int> Al(27,13);
   std::pair<int,int> Fe(56,26);
   std::pair<int,int> Pb(207,82);
   
   std::pair<int,int> t(3,1);
   std::pair<int,int> alpha(4,2);
   std::pair<int,int> Li6(6,3);
   std::pair<int,int> Li7(7,3);
   std::pair<int,int> C(12,6);
   

   materials.push_back(H);
   materials.push_back(d);
   materials.push_back(O);
   materials.push_back(Al);
   materials.push_back(Fe);
   materials.push_back(Pb);

   projectiles.push_back(d);
   projectiles.push_back(t);
   projectiles.push_back(alpha);
   projectiles.push_back(Li6);
   projectiles.push_back(Li7);
   projectiles.push_back(C);
   
   energies.push_back(100.);
   energies.push_back(250.);
   energies.push_back(500.);
   energies.push_back(800.);

   std::vector<std::pair<int,int> >::iterator  imat;
   std::vector<std::pair<int,int> >::iterator ipro;
   std::vector<double>::iterator iener;

   for (imat=materials.begin();imat<materials.end();imat++)
   {
       for (ipro=projectiles.begin();ipro<projectiles.end();ipro++)
       {
           for (iener=energies.begin();iener<energies.end();iener++)
	   {
	        int m_a=imat->first;
		int m_z=imat->second;
		int p_a=ipro->first;
		int p_z=ipro->second;
		double e=*iener;
		std::strstream f_name;
		f_name << m_a <<"_" << m_z<<"_" << p_a<<"_"<< p_z<<
		"_"<<e<<".in";
		std::cout << m_a <<" " << m_z<<" " << e<< "f_name " <<
		f_name.str() << std::endl;
	   	std::ofstream  file(f_name.str());
	   	file << m_a <<" " << m_z<<" " << std::endl;
		file << "1" << std::endl;
		file << p_a <<" " << p_z<<" " << std::endl;
		file << "20000" << std::endl;   // number of events
		file << e << std::endl;
		file << "1" << std::endl;
		file.close();
		
	   }
	   
       }
   }
}   
