// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4MEMSTAT_HH_
#define G4MEMSTAT_HH_

#include "globals.hh"

namespace G4MemStat
{

	struct MemStat;

	std::ostream & operator<<(std::ostream &os, const MemStat& p);

	struct MemStat
	{

        friend std::ostream & operator<<(std::ostream &os, const MemStat& p);
		double vmz;
		double mem;

		MemStat() : vmz(0), mem(0)
		{;}
		MemStat(const MemStat& right)
		{
			vmz = right.vmz;
			mem = right.mem;
		}

		MemStat operator-(const MemStat& right)
		{
			MemStat output;
			output.vmz = this->vmz-right.vmz;
			output.mem = this->mem-right.mem;
			return output;
		}
	};

	MemStat MemoryUsage();

	std::ostream & operator<<(std::ostream &os, const MemStat& p);
}
#endif /* MEMSTAT_HH_ */
