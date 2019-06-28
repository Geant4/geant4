// -------------------------------------------------------------------
// $Id: plot.C 70323 2013-05-29 07:57:44Z gcosmo $
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// *********************************************************************

#include <iostream>
#include "TROOT.h"
#include "../dnadamage1/include/DNAVolumeType.hh"
using namespace std;

template<typename T>  
class ThreeVector
{
private:
    T _x, _y, _z;
public:
    ThreeVector():_x(0),_y(0),_z(0){}
    ThreeVector(T x, T y, T z)
        :_x(x),_y(y),_z(z){}
    ~ThreeVector(){}
    T x() const
    {
        return _x;
    }
    T y() const
    {
        return _y;
    }
    T z() const
    {
    return _z;
    }
    
    bool operator ==(const ThreeVector<T>& right) const
    {
        return (_x == right._x) && 
               (_y == right._y) &&
               (_z == right._z);
    }
            
    ThreeVector<T>& operator =(const ThreeVector<T>& right) = default;

ClassDef(ThreeVector,1)
};

#if !defined(__CLING__)
ClassImp(ThreeVector);
#endif

class Molecule
{
public:
    Molecule(){}
    Molecule(string name, 
             int copyNumber, 
             const ThreeVector<double>& position, 
             int strand)
             : fName(name)
             , fCopyNumber(copyNumber)
             , fPosition(position)
             , fStrand(strand)
            {}
    ~Molecule(){}
public:
    string fName;
    string fMaterial;
    int fCopyNumber;
    int fStrand;

    ThreeVector<double> fPosition;

    double fRadius;
    double fRadiusWater;

    ClassDef(Molecule,1)
};

#if !defined(__CLING__)
ClassImp(Molecule);
#endif

std::vector<Molecule> molecule()
{
    std::vector<Molecule> fMolecules;
    double size;
    string name;
    ifstream file("VoxelStraight.fab2g4dna");
    if(!file.is_open())
    {
        string msg ="VoxelStraight.fab2g4dna could not be opened";
        throw std::invalid_argument(msg);
    }

   string line;
    while(getline(file, line) )
    {
        if(line.empty()) 
        {
            continue;
        }
         
        istringstream issLine(line);
        string firstItem;
        issLine >> firstItem;
        if("_Size" == firstItem)
        {
            issLine >> size;
        }
        else if("_pl" == firstItem)
        {
            string name;
            issLine >> name;

            string material;
            issLine >> material;

            int strand;
            issLine >> strand;

            int copyNumber;
            issLine >> copyNumber;

            double x;
            issLine >> x;

            double y;
            issLine >> y;

            double z;
            issLine >> z;

            Molecule molecule(name, 
                              copyNumber, 
                              ThreeVector<double>(x, y, z), 
                              strand);
            fMolecules.push_back(molecule);
        }
    }
    file.close();
    
    return fMolecules;
}
