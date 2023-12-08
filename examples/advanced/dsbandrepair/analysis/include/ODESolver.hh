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
/// \file ODESolver.hh
/// \brief Definition of the ODESolver class

#ifndef ODESolver_h
#define ODESolver_h 1

#include <vector>
#include <map>
#include <functional>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

using myODEs = std::vector<double>(*)(double ,std::vector<double>) ;
std::vector<double> operator*(const std::vector<double> v, double alfa);
std::vector<double> operator+(const std::vector<double> v, double alfa);
std::vector<double> operator+(const std::vector<double> v1, const std::vector<double> v2);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class ODESolver
{
public:
    ODESolver();
    ~ODESolver() = default;
    void Embedded_RungeKutta_Fehlberg(  
        std::function<std::vector<double>(double,std::vector<double>)>, std::vector<double> &y,
        double start,double end,double stepsize=-1, double epsilon = 1e-3,
        std::vector<double> *time_observer=nullptr,std::vector<std::vector<double>> *state_observer=nullptr);
    void SetNstepsForObserver(unsigned int ndt) {fNstepsForObserver = ndt;}
    void RungeKutta4(  
        std::function<std::vector<double>(double,std::vector<double>)>, std::vector<double> &y,
        double start,double end,double stepsize=-1, 
        std::vector<double> *time_observer=nullptr,std::vector<std::vector<double>> *state_observer=nullptr);
private:
    double RungeKutta_Fehlberg(std::function<std::vector<double>(double,std::vector<double>)>, 
            std::vector<double> &y,double t, double stepsize=-1);
    void absValuesVector(std::vector<double> &vIn);
    unsigned int fNstepsForObserver{1};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

inline void ODESolver::absValuesVector(std::vector<double> &vIn)
{
    for (double &val : vIn) {
        if (val < 0) val *= -1.;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif