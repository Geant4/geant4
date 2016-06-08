
#include "G4Solver.hh"

template <class Function> 
G4Solver<Function>::G4Solver(const G4Solver & right)
{
	MaxIter = right.MaxIter;
	tolerance = right.tolerance;
	a = right.a;
	b = right.b;
	root = right.root;
}

// operators
template <class Function> 
G4Solver<Function> & G4Solver<Function>::operator=(const G4Solver & right)
{
	MaxIter = right.MaxIter;
	tolerance = right.tolerance;
	a = right.a;
	b = right.b;
	root = right.root;	
	return *this;
}

template <class Function> 
G4bool G4Solver<Function>::operator==(const G4Solver & right) const
{
	if (this == &right) return true;
	else return false;
}

template <class Function> 
G4bool G4Solver<Function>::operator!=(const G4Solver & right) const
{
	return !operator==(right);
}
	



