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
/*
 * G4Serialize.cc
 *
 *  Created on: Jul 10, 2015
 *      Author: mkaramit
 */


#include "G4Serialize.hh"

template<>
void WRITE<G4String>(std::ostream& out, const G4String& name)
{
  size_t size = name.size();
  out.write((char*)(&size), sizeof(size_t));
  out.write(name.c_str(), size);
}

//_____________________________________________________________________________

class UReadBinaryString
{
public:
    static G4String read(std::istream &is, size_t size)
    {
      G4String  returnStr;
        if(size > 0)
        {
            char* buff = new char[size];
            is.read(buff, size);
            returnStr.assign(buff, size);
            delete [] buff;
        }

        return returnStr;
    }
};

//_____________________________________________________________________________

template<>
void READ<G4String>(std::istream& in, G4String& name)
{
  size_t size;
  in.read((char*)(&size), sizeof(size_t));
  name = UReadBinaryString::read(in, size);
}
