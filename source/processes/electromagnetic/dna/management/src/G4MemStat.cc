/*
 * MemStat.cc
 *
 *  Created on: 4 févr. 2014
 *      Author: kara
 */

#include "G4MemStat.hh"

#if ( defined(__MACH__) && defined(__clang__) && defined(__x86_64__) ) || \
    ( defined(__MACH__) && defined(__GNUC__) && __GNUC__>=4 && __GNUC_MINOR__>=7 ) || \
    defined(__linux__) || defined(_AIX)

#include <unistd.h>

#endif

#include <ios>
#include <iostream>
#include <fstream>
#include <string>

namespace G4MemStat
{

using std::ios_base;
using std::ifstream;
using std::string;

MemStat MemoryUsage()
{
	MemStat output;

#if ( defined(__MACH__) && defined(__clang__) && defined(__x86_64__) ) || \
    ( defined(__MACH__) && defined(__GNUC__) && __GNUC__>=4 && __GNUC_MINOR__>=7 ) || \
    defined(__linux__) || defined(_AIX)

	// 'file' stat seems to give the most reliable results
	//
	ifstream stat_stream("/proc/self/stat",ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	//
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	>> utime >> stime >> cutime >> cstime >> priority >> nice
	>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	stat_stream.close();

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	output.vmz = vsize / 1024.0;
	output.mem = rss * page_size_kb;
#endif

	return output;
}

std::ostream & operator<<(std::ostream &os, const MemStat& memStat)
{
    return os << "( vmz: "<< memStat.vmz << ", " << "mem: " << memStat.mem << ")" ;
}

}
