#!/usr/bin/env python3

"""Check Geant4's module dependency graph for cycles

Geant4 organises source code into "Modules", with each Module being a
directory containing headers and sources in `include/` and `src/`
subdirectories respectively. One or more Modules are grouped/compiled
into the actual libraries.

CMake will raise errors if there are circular dependencies between
shared libraries (cycles are allowed between static libraries). However,
CMake cannot detect cycles between Modules directly, are these are a
Geant4 construct. Whilst cycles between modules that get added to the
same library are technically o.k., they still indicate either:

- A bad design/organisation of code
- A mistake in the declared dependencies for a module(s)

To help developers and CI pick module cycles, Geant4's CMake scripts
write out the declared module dependencies in a text file as an Adjacency
List graph. This comprises one line per module, with the first column
being the module name, and any remaining columns being the modules the
first one "depends on", i.e. uses code from.

This program reads that file and uses Python's graphlib module to try
topologically sorting this graph. Failure to sort indicates a cycle, and
this will be printed to stdout and a non-zero return code emitted on exit.
"""

import argparse
import graphlib
import sys

if __name__ == "__main__":
    # Parse command line to get file to load
    parser = argparse.ArgumentParser(
        description=str(__doc__), formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-f",
        "--file",
        required=True,
        metavar="FILE",
        help="file to read adjacency list from",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="print topologically sorted list of modules",
    )
    args = parser.parse_args()

    try:
        # 1. Try and construct graph from input adjancey list file
        adjList = {}
        with open(args.file) as f:
            for line in f:
                if not line.strip().startswith("#"):
                    nodes = line.rstrip("\n").split(" ")
                    adjList[nodes[0]] = nodes[1:]

        if not adjList:
            print(f"warning: graph read from '{args.file}' is empty")

        # NB also possible with networkx, but requires pip install:
        # G = nx.read_adjlist(args.file, create_using=nx.DiGraph)
        # cycles = nx.find_cycle(G, orientation='original')

        # Topo sort throws cycle error if one occurs during prepare
        ts = graphlib.TopologicalSorter(adjList)
        ts.prepare()

        # Verbose print of nodes in topological order
        # NB: This uses the full graphlib interface in case we want to
        # print out nodes grouped by level in the graph.
        if args.verbose:
            while ts.is_active():
                nodes = ts.get_ready()
                for n in nodes:
                    print(n)
                ts.done(*nodes)

        print(f"pass: No cycles detected in graph defined in '{args.file}'")

    except OSError as err:
        print(f"OS error: {err}")
        sys.exit(1)
    except graphlib.CycleError as err:
        print(
            f"error: cycles detected in input graph from '{args.file}':",
            file=sys.stderr,
        )
        cycle = " <- ".join(err.args[1])
        print(f"cycle: {cycle}", file=sys.stderr)
        sys.exit(1)
