#!/usr/bin/env python3

"""Check Geant4 source code modules

Geant4 organises source code into "Modules", with each Module being a
directory containing headers and sources in `include/` and `src/`
subdirectories respectively. One or more Modules are grouped/compiled
into the actual libraries.

This program can be used to query Modules as declared to the build system to:

- List all declared modules
- Print the list of public headers provided by a module
- Print modules required by a module based on inclusion of headers
- Check consistency between declared dependencies of a module and
  those its sources actually use via inclusion of headers
- Check for cycles in the graph of declared module dependencies
"""

import argparse
import csv
import os
import sys
import re
import json
import graphlib


def initdb(filename):
    db = {}
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            db[row[0]] = {'location': row[1],
                          'headers': set(row[2].split(';')) - set(['']),
                          'public_deps': set(row[3].split(';')) - set(['']),
                          'private_deps': set(row[4].split(';')) - set(['']),
                          'interface_deps': set(row[5].split(';')) - set([''])}

    return db


def scan_includes(input_file):
    result = set()
    lines = open(input_file).readlines()
    for line in lines:
        m = re.search("^ *# *include *\"(.*)\"", line)
        if m:
            result.add(m.group(1))

    return result


def scan_files(input_path):
    used_headers = set()
    with os.scandir(input_path) as hdr_it:
        for entry in hdr_it:
            used_headers |= scan_includes(entry.path)

    return used_headers


def find_modules(list_of_headers, module_db):
    found_modules = set()
    orphan_hdrs = set()
    for h in list_of_headers:
        ms = [k for k, v in module_db.items() if h in v['headers']]
        found_modules |= set(ms)
        if not ms:
            orphan_hdrs.add(h)

    return {'modules': found_modules, 'headers': orphan_hdrs}


def usage_requirements(module_name, module_db):
    """ Find modules needed publically and privately by input module
    """
    module_path = module_db[module_name]['location']
    module_headers = module_db[module_name]['headers']

    # Determine public usage reqs
    includes_from_headers = scan_files(os.path.join(module_path, "include"))
    includes_from_headers -= module_headers
    public_deps = find_modules(includes_from_headers, module_db)

    # The same for private
    includes_from_srcs = scan_files(os.path.join(module_path, "src"))
    includes_from_srcs -= module_headers
    private_deps = find_modules(includes_from_srcs, module_db)
    # Public deps are higher priority
    private_deps['modules'] -= public_deps['modules']
    private_deps['headers'] -= public_deps['headers']

    # Transform results to output dict
    d = {'module': module_name,
         'dependencies': {
             'public': sorted(public_deps['modules']),
             'private': sorted(private_deps['modules'])
         },
         'external_headers': {
             'public': sorted(public_deps['headers']),
             'private': sorted(private_deps['headers'])
         }
         }
    return d


def check_consistency(module_name, module_db):
    """ Check module declared/apparent dependencies for consistency
    """
    # NB: Can have false positives from externals G4expat, G4clhep, G4zlib, G4tools (plus imported :: targets)
    def filter_dependencies(dep_list):
        return set([x for x in dep_list if not re.match("^.+::.+", x)]) - set(['G4expat', 'G4clhep', 'G4zlib', 'G4tools', 'G4ptl'])

    ur = usage_requirements(module_name, module_db)
    apparent_public_deps = filter_dependencies(ur['dependencies']['public'])
    apparent_private_deps = filter_dependencies(ur['dependencies']['private'])

    declared_public_deps = filter_dependencies(
        module_db[module_name]['public_deps'])
    declared_private_deps = filter_dependencies(
        module_db[module_name]['private_deps'])

    # Collate any consistency errors
    report = []
    # - Declared dependencies duplicated between public/private
    #   Later checks will help distinguish what to do with this
    duplicated_deps = declared_public_deps & declared_private_deps
    if duplicated_deps:
        report.append(
            f'- has duplicated PUBLIC/PRIVATE dependencies: {duplicated_deps}')

    # - Apparent dep not declared
    missing_public_deps = apparent_public_deps - declared_public_deps
    if missing_public_deps:
        report.append(
            f'- may require PUBLIC dependencies: {missing_public_deps}')

    missing_private_deps = apparent_private_deps - declared_private_deps
    if missing_private_deps:
        report.append(
            f'- may require PRIVATE dependencies: {missing_private_deps}')

    # - Declared dep not in apparent (overdeclared)
    overdeclared_public_deps = declared_public_deps - apparent_public_deps
    if overdeclared_public_deps:
        report.append(
            f'- may not require PUBLIC dependencies: {overdeclared_public_deps}')

    overdeclared_private_deps = declared_private_deps - apparent_private_deps
    if overdeclared_private_deps:
        report.append(
            f'- may not require PRIVATE dependencies: {overdeclared_private_deps}')

    return report


def find_cycles(module_db, verbose):
    """ Check for any cycles in the complete module dependency graph
    """
    # Build adjacency list
    adjlist = {}
    for k, v in module_db.items():
        adjlist[k] = list(v['public_deps'] | v['private_deps']
                          | v['interface_deps'])

    try:
        # Topo sort throws cycle error if one occurs during prepare
        ts = graphlib.TopologicalSorter(adjlist)
        ts.prepare()

        # Verbose print of nodes in topological order
        # NB: This uses the full graphlib interface in case we want to
        # print out nodes grouped by level in the graph.
        if verbose:
            generation = 0
            while ts.is_active():
                nodes = ts.get_ready()
                print(f'Generation: {generation}')
                for n in nodes:
                    print(f'  {n}')
                ts.done(*nodes)
                generation += 1

        print(f"pass: No cycles detected in module dependency graph")
    except graphlib.CycleError as err:
        print(
            f"error: cycles detected in module dependency graph",
            file=sys.stderr,
        )
        cycle = " <- ".join(err.args[1])
        print(f"cycle: {cycle}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    # Parse command line to get command to run
    parser = argparse.ArgumentParser(
        description=str(__doc__), formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-db", default='G4ModuleInterfaceMap.csv',
                        metavar="FILE", help="module interface map file")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="verbose output")

    query_group = parser.add_mutually_exclusive_group(required=True)

    query_group.add_argument(
        "-l", "--list",
        action='store_true',
        help="list declared source code modules"
    )
    query_group.add_argument(
        "-i", "--interface",
        metavar="<module>",
        help="print public headers of module",
    )
    query_group.add_argument(
        "-p", "--provides",
        metavar="<header>",
        help="print module that provides this header"
    )
    query_group.add_argument(
        "-u", "--usage-requirements",
        metavar="<module>",
        help="determine and print usage requirements of module"
    )
    query_group.add_argument(
        "-c", "--check-consistency",
        metavar="<module>",
        help="check declared and apparent module dependencies for basic consistency"
    )
    query_group.add_argument(
        "--find-cycles",
        action="store_true",
        help="find cycles in graph of declared modules dependencies"
    )
    query_group.add_argument(
        "--find-inconsistencies",
        action="store_true",
        help="find inconsistencies in apparent/declared dependencies of all modules"
    )
    args = parser.parse_args()

    # Initialize the module/header "database"
    try:
        db = initdb(args.db)
    except OSError as err:
        print(f"Could not initalize module DB: {err}", file=sys.stderr)
        sys.exit(1)

    # Implementations
    if args.list:
        print("\n".join(db.keys()))
    elif args.interface:
        try:
            print("\n".join(db[args.interface]['headers']))
        except KeyError as err:
            print(f"No module named {err}", file=sys.stderr)
            sys.exit(1)
    elif args.provides:
        # find the header, there should not be duplicates!
        mods = [k for k, v in db.items() if args.provides in v['headers']]
        if len(mods) == 0:
            print(
                f"No module provides header '{args.provides}'",
                file=sys.stderr)
            sys.exit(1)
        if len(mods) > 1:
            print(
                f"Header '{args.provides}' is provided by multiple modules: {mods}",
                file=sys.stderr)
            sys.exit(1)

        print(mods[0])
    elif args.usage_requirements:
        ur = usage_requirements(args.usage_requirements, db)
        print(json.dumps(ur, indent=2))
    elif args.check_consistency:
        cc = check_consistency(args.check_consistency, db)
        if cc:
            print(f'{args.check_consistency}:', file=sys.stderr)
            print("\n".join(cc), file=sys.stderr)
            sys.exit(1)
        else:
            print(
                f'{args.check_consistency} appears consistent')
    elif args.find_cycles:
        find_cycles(db, args.verbose)
    elif args.find_inconsistencies:
        inconsistent = {}
        for m in db.keys():
            cc = check_consistency(m, db)
            if cc:
                inconsistent[m] = "\n".join(cc)

        if len(inconsistent) > 0:
            print(f'Inconsistent dependencies in modules:', file=sys.stderr)
            for k, v in inconsistent.items():
                print(f'{k}:\n{v}', file=sys.stderr)

            sys.exit(1)
