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
/// \file ODESolver.cc
/// \brief Implementation of the ODESolver class

#include "ODESolver.hh"
#include <iostream>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector<double> operator*(const std::vector<double> v, double alfa)
{
	std::vector<double> vout;
	for (auto const val : v) vout.push_back(val*alfa);
	return vout;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector<double> operator+(const std::vector<double> v, double alfa)
{
	std::vector<double> vout;
	for (auto const val : v) vout.push_back(val + alfa);
	return vout;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector<double> operator+(const std::vector<double> v1, const std::vector<double> v2)
{
	std::vector<double> vout;
	for (size_t i=0;i<v1.size();i++) vout.push_back(v1.at(i) + v2.at(i));
	return vout;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ODESolver::ODESolver(): fNstepsForObserver(1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

double ODESolver::RungeKutta_Fehlberg( std::function<std::vector<double>(double,std::vector<double>)> 
func,std::vector<double> &y, double t, double stepsize)
{
	//based on https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
	const int nk=6;
	double h = stepsize;
	double CH[nk]={47./450.,0,12./25.,32./255.,1./30.,6./25.};
	double CT[nk]={-1./150.,0.,3./100.,-16./75.,-1./20.,6./25.};
	double A[nk]={0.,2./9.,1./3.,3./4.,1.,5./6.};
	double B21=2./9., B31=1./12., B41=69./128., B51=-17./12., B61=65./432.;
	double B32=1./4., B42=-243./128., B52=27./5., B62=13./16.;
	double B43=135./64., B53=-27./5., B63=13./16.;
	double B54=16./15., B64=4./27.;
	double B65=5./144.;
	double maxError = 1.;
	std::vector<double> k1 = func(t+A[0]*h,y)*h;
	std::vector<double> k2 = func(t+A[1]*h,y + k1*B21)*h;
	std::vector<double> k3 = func(t+A[2]*h,y + k1*B31 + k2*B32)*h;
	std::vector<double> k4 = func(t+A[3]*h,y + k1*B41 + k2*B42 + k3*B43)*h;
	std::vector<double> k5 = func(t+A[4]*h,y + k1*B51 + k2*B52 + k3*B53 + k4*B54)*h;
	std::vector<double> k6 = func(t+A[5]*h,y + k1*B61 + k2*B62 + k3*B63 + k4*B64 + k5*B65)*h;
	y = y + k1*CH[0] + k2*CH[1] + k3*CH[2] + k4*CH[3] + k5*CH[4] + k6*CH[5];
	auto TE = k1*CT[0] + k2*CT[1] + k3*CT[2] + k4*CT[3] + k5*CT[4] + k6*CT[5];
	absValuesVector(TE);
	maxError = *std::max_element(TE.begin(),TE.end());
	
	k1.clear(); k1.shrink_to_fit();
	k2.clear(); k2.shrink_to_fit();
	k3.clear(); k3.shrink_to_fit();
	k4.clear(); k4.shrink_to_fit();
	k5.clear(); k5.shrink_to_fit();
	k6.clear(); k6.shrink_to_fit();
	TE.clear(); TE.shrink_to_fit();
	return  maxError;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ODESolver::Embedded_RungeKutta_Fehlberg(	
	std::function<std::vector<double>(double,std::vector<double>)> func, std::vector<double> &y,
	double start,double end,double stepsize,double epsilon,
	std::vector<double> *time_observer,std::vector<std::vector<double>> *state_observer)
{
	double t = start;
	double h = stepsize;
	int nsteps = 0;
	if (h < 0) h = (end -  start)/(10000.);	
	if (time_observer) time_observer->push_back(t);
	if (state_observer) state_observer->push_back(y);
	auto ytemp = y;
	while (t < end)
	{
		ytemp = y;
		double maxerror = RungeKutta_Fehlberg(func,ytemp,t,h);
		double scale = 0.9*std::pow(epsilon/maxerror,1./5.);
		
		double hnew = h*scale;
		while (maxerror > epsilon)
		{
			ytemp = y;
			maxerror = RungeKutta_Fehlberg(func,ytemp,t,hnew);
			scale = 0.9*std::pow(epsilon/maxerror,1./5.);
			hnew = hnew*scale;
		}
		h = hnew;
		y = ytemp;
		t += h;
		if (t > end) break;
		if ( time_observer || state_observer) {
			nsteps++;
			if (nsteps%fNstepsForObserver == 0) {
				if (time_observer) time_observer->push_back(t);
				if (state_observer) state_observer->push_back(y);	
				nsteps = 0;
			}
		}
	}

	ytemp.clear(); ytemp.shrink_to_fit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ODESolver::RungeKutta4(  
	std::function<std::vector<double>(double,std::vector<double>)> func, std::vector<double> &y,
    double start,double end,double stepsize, 
    std::vector<double> *time_observer,std::vector<std::vector<double>> *state_observer)
{
	double t = start;
	double h = stepsize;
	int nsteps = 0;
	if (h < 0) h = (end -  start)/(10000.);	
	if (time_observer) time_observer->push_back(t);
	if (state_observer) state_observer->push_back(y);
	while (t < end)
	{
		std::vector<double> k1 = func(t,y)*h;
		std::vector<double> k2 = func(t+0.5*h,y + k1*0.5)*h;
		std::vector<double> k3 = func(t+0.5*h,y + k2*0.5)*h;
		std::vector<double> k4 = func(t+h,y + k3)*h;
		t += h;
		if (t > end) break;
		y = y +(k1 +k2*2+k3*2+k4)*(1./6.0);
		if ( time_observer || state_observer) {
			nsteps++;
			if (nsteps%fNstepsForObserver == 0) {
				if (time_observer) time_observer->push_back(t);
				if (state_observer) state_observer->push_back(y);	
				nsteps = 0;
			}
		}
		k1.clear(); k1.shrink_to_fit();
		k2.clear(); k2.shrink_to_fit();
		k3.clear(); k3.shrink_to_fit();
		k4.clear(); k4.shrink_to_fit();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....