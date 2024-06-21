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
from collections import defaultdict, OrderedDict


def initdb(filename):
    db = {'modules': {}, 'build_settings': set()}
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            db['modules'][row[0]] = {'location': row[1],
                                     'headers': set(row[2].split(';')) - set(['']),
                                     'private_headers': set(row[3].split(';')) - set(['']),
                                     'sources': set(row[4].split(';')) - set(['']),
                                     'public_deps': set(row[5].split(';')) - set(['']),
                                     'private_deps': set(row[6].split(';')) - set(['']),
                                     'interface_deps': set(row[7].split(';')) - set(['']),
                                     'parent_target': row[8]}
            db['build_settings'].add(row[9])

    return db


def get_modules(db):
    return db['modules'].keys()


def get_module(name, db):
    return db['modules'][name]


def get_location(name, db):
    return get_module(name, db)['location']


def get_header_path(name, db):
    return os.path.join(get_location(name, db), 'include')


def get_private_header_path(name, db):
    return os.path.join(get_location(name, db), 'include/private')


def get_source_path(name, db):
    return os.path.join(get_location(name, db), 'src')


def get_headers(name, db):
    return get_module(name, db)['headers']


def get_headers_absolute(name, db):
    base_path = get_header_path(name, db)
    return set(os.path.join(base_path, h) for h in get_headers(name, db))


def get_private_headers(name, db):
    return get_module(name, db)['private_headers']


def get_private_headers_absolute(name, db):
    base_path = get_private_header_path(name, db)
    return set(os.path.join(base_path, h) for h in get_private_headers(name, db))


def get_sources(name, db):
    return get_module(name, db)['sources']


def get_sources_absolute(name, db):
    base_path = get_source_path(name, db)
    return set(os.path.join(base_path, h) for h in get_sources(name, db))


def get_public_deps(name, db):
    return get_module(name, db)['public_deps']


def get_private_deps(name, db):
    return get_module(name, db)['private_deps']


def get_interface_deps(name, db):
    return get_module(name, db)['interface_deps']


def get_parent_target(name, db):
    return get_module(name, db)['parent_target']


def what_provides(h, db):
    return set([m for m in get_modules(db) if h in get_headers(m, db)])


def what_provides_private(h, db):
    return set([m for m in get_modules(db) if h in get_private_headers(m, db)])


def what_requires(m, db):
    return set([om for om in get_modules(db) if m in get_public_deps(om, db)])


def what_requires_private(m, db):
    return set([om for om in get_modules(db) if m in get_private_deps(om, db)])


def has_setting(setting, build_settings):
    if setting.startswith('!'):
        return setting[1:] not in build_settings
    else:
        return setting in build_settings


def scan_includes(input_file, settings):
    result = set()
    lines = open(input_file).readlines()
    for line in lines:
        m = re.search("^ *# *include *\"(.*)\"(.*)$", line)
        if m:
            c = re.search(
                "\\s*//\\s*no_geant4_module_check\\s*(\\((.*)\\))?.*$", m.group(2))
            if c:
                conds = set()
                if c.group(2):
                    conds = set(
                        filter(None, [''.join(c.split(' ')) for c in c.group(2).strip().split(',')]))

                if not conds or any(has_setting(s, settings) for s in conds):
                    continue

            result.add(m.group(1))

    return result


def scan_files(files, settings):
    used_headers = set()
    for f in files:
        used_headers |= scan_includes(f, settings)

    return used_headers


def scan_module_headers(name, db, full_scan):
    if full_scan:
        # NB: know that will have private/ subdir, so only consider files
        header_list = [e.path for e in os.scandir(
            get_header_path(name, db)) if e.is_file()]
    else:
        header_list = get_headers_absolute(name, db)

    return scan_files(header_list, db['build_settings'])


def scan_module_private_headers(name, db, full_scan):
    if full_scan:
        header_list = [e.path for e in os.scandir(
            get_private_header_path(name, db))]
    else:
        header_list = get_private_headers_absolute(name, db)

    return scan_files(header_list, db['build_settings'])


def scan_module_sources(name, db, full_scan):
    if full_scan:
        source_list = [e.path for e in os.scandir(get_source_path(name, db))]
    else:
        source_list = get_sources_absolute(name, db)

    return scan_files(source_list, db['build_settings'])


def find_modules(list_of_headers, module_db):
    found_modules = set()
    blame_headers = defaultdict(set)
    orphan_hdrs = set()
    for h in list_of_headers:
        ms = what_provides(h, module_db)
        found_modules |= ms
        if not ms:
            orphan_hdrs.add(h)
        else:
            blame_headers[ms.pop()].add(h)

    return {'modules': found_modules, 'blame': {k: sorted(v) for k, v in blame_headers.items()}, 'headers': orphan_hdrs}


def find_modules_private(list_of_headers, module_db):
    found_modules = set()
    blame_headers = defaultdict(set)
    orphan_hdrs = set()
    for h in list_of_headers:
        ms = what_provides_private(h, module_db)
        found_modules |= ms
        if not ms:
            orphan_hdrs.add(h)
        else:
            blame_headers[ms.pop()].add(h)

    return {'modules': found_modules, 'blame': {k: sorted(v) for k, v in blame_headers.items()}, 'headers': orphan_hdrs}


def usage_requirements(module_name, module_db, full_scan):
    """ Find modules needed publically and privately by input module
    """
    try:
        module_path = get_location(module_name, module_db)
        module_headers = get_headers(module_name, module_db)
        module_private_headers = get_private_headers(module_name, module_db)
    except KeyError as err:
        print(f'No module named \'{module_name}\'', file=sys.stderr)
        sys.exit(1)

    # Determine public usage reqs
    includes_from_headers = scan_module_headers(
        module_name, module_db, full_scan)
    includes_from_headers -= module_headers
    public_deps = find_modules(includes_from_headers, module_db)
    access_violation_public = find_modules_private(
        includes_from_headers, module_db)

    # Determine private usage reqs
    # private headers...
    try:
        includes_from_private_headers = scan_module_private_headers(
            module_name, module_db, full_scan)
        includes_from_private_headers -= module_headers
        includes_from_private_headers -= module_private_headers
    except FileNotFoundError:
        includes_from_private_headers = set()

    # ... and sources if they exist
    try:
        includes_from_srcs = scan_module_sources(
            module_name, module_db, full_scan)
        includes_from_srcs -= module_headers
        includes_from_srcs -= module_private_headers
    except FileNotFoundError:
        includes_from_srcs = set()

    # Join private/src and remove public deps, which are higher priority
    # We don't do this for blames as these can be unique
    private_deps = find_modules(
        includes_from_private_headers | includes_from_srcs, module_db)
    private_deps['modules'] -= public_deps['modules']
    private_deps['headers'] -= public_deps['headers']

    # Determine any access violations
    access_violation_private = find_modules_private(
        includes_from_private_headers | includes_from_srcs, module_db)

    # Transform results to output dict
    d = {'module': module_name,
         'dependencies': {
             'public': sorted(public_deps['modules']),
             'private': sorted(private_deps['modules'])
         },
         'blame': {
             'public': OrderedDict(sorted(public_deps['blame'].items())),
             'private': OrderedDict(sorted(private_deps['blame'].items()))
         },
         'external_headers': {
             'public': sorted(public_deps['headers']),
             'private': sorted(private_deps['headers'])
         },
         'access_violations': {
             'public': sorted(access_violation_public['modules']),
             'private': sorted(access_violation_private['modules'])
         }
         }
    return d


def check_consistency(module_name, module_db, full_scan):
    """ Check module declared/apparent dependencies for consistency
    """
    # NB: Can have false positives from externals G4expat, G4clhep, G4zlib, G4tools (plus imported :: targets and stdlib)
    def filter_dependencies(dep_list, module_name):
        return set([x for x in dep_list if not re.match("^.+::.+", x)]) - set(['G4expat', 'G4clhep', 'G4zlib', 'G4tools', 'G4ptl', 'stdc++fs'])

    ur = usage_requirements(module_name, module_db, full_scan)
    apparent_public_deps = filter_dependencies(
        ur['dependencies']['public'], module_name)
    apparent_private_deps = filter_dependencies(
        ur['dependencies']['private'], module_name)

    declared_public_deps = filter_dependencies(
        get_public_deps(module_name, module_db), module_name)
    declared_interface_deps = filter_dependencies(
        get_interface_deps(module_name, module_db), module_name)
    declared_private_deps = filter_dependencies(
        get_private_deps(module_name, module_db), module_name)

    public_access_violations = filter_dependencies(
        ur['access_violations']['public'], module_name)
    private_access_violations = filter_dependencies(
        ur['access_violations']['private'], module_name)

    # Collate any consistency errors
    report = []
    # - Declared dependencies duplicated between public/private
    #   Later checks will help distinguish what to do with this
    duplicated_deps = declared_public_deps & declared_private_deps
    if duplicated_deps:
        report.append(
            f'  - has duplicated PUBLIC/PRIVATE dependencies: {duplicated_deps}')

    # - Apparent dep not declared
    missing_public_deps = apparent_public_deps - \
        declared_public_deps - declared_interface_deps
    if missing_public_deps:
        report.append(
            f'  + may require PUBLIC or INTERFACE dependencies: {missing_public_deps}')

    missing_private_deps = apparent_private_deps - declared_private_deps
    if missing_private_deps:
        report.append(
            f'  + may require PRIVATE dependencies: {missing_private_deps}')

    # - Declared dep not in apparent (overdeclared)
    overdeclared_public_deps = declared_public_deps - apparent_public_deps
    if overdeclared_public_deps:
        report.append(
            f'  - may not require PUBLIC dependencies: {overdeclared_public_deps}')

    overdeclared_private_deps = declared_private_deps - apparent_private_deps
    if overdeclared_private_deps:
        report.append(
            f'  - may not require PRIVATE dependencies: {overdeclared_private_deps}')

    # - Module exposes a header that's declared private locally or in another module
    if public_access_violations:
        report.append(
            f'  ! public interface #include-s PRIVATE headers from modules: {public_access_violations}')

    # - Module uses a header that's declared private
    if private_access_violations:
        report.append(
            f'  ! implementation uses PRIVATE headers from modules: {private_access_violations}')

    return report


def do_provides(header_name, module_db, verbose):
    # find the header, there should not be duplicates!
    mods = what_provides(header_name, module_db)
    if len(mods) == 0:
        print(f"No module provides header '{header_name}'", file=sys.stderr)
        sys.exit(1)
    if len(mods) > 1:
        print(
            f"Header '{args.provides}' is provided by multiple modules: {mods}", file=sys.stderr)
        sys.exit(1)

    print(mods.pop())


def do_requires(module_name, module_db):
    """find all modules that use module_name
    """
    try:
        public_users = sorted(what_requires(module_name, module_db))
        private_users = sorted(what_requires_private(module_name, module_db))
    except KeyError as err:
        print(f'No module named \'{module_name}\'', file=sys.stderr)
        sys.exit(1)

    print(f'{module_name} is required by:')

    if len(public_users) != 0:
        print('- PUBLIC:')
        print("  -", "\n  - ".join(public_users))

    if len(private_users) != 0:
        print('PRIVATE:')
        print("  -", "\n  - ".join(private_users))


def do_count_usage(module_db, verbose):
    """print modules ordered by number of modules linking to them"""
    reqcount = dict.fromkeys(get_modules(module_db), 0)

    for m in get_modules(module_db):
        reqcount[m] += len(what_requires(m, module_db))
        reqcount[m] += len(what_requires_private(m, module_db))

    for m in sorted(reqcount, key=reqcount.get, reverse=True):
        if verbose:
            print(m, f'({reqcount[m]} users)')
        else:
            print(m)


def do_check_consistency(module_name, module_db, full_scan):
    try:
        cc = check_consistency(module_name, module_db, full_scan)
        if cc:
            print(f'{module_name} has inconsistent dependencies:', file=sys.stderr)
            print("\n".join(cc), file=sys.stderr)
            sys.exit(1)
        else:
            print(f'Module {module_name} appears consistent')
    except KeyError as err:
        print(f'No module named \'{module_name}\'', file=sys.stderr)
        sys.exit(1)


def do_find_cycles(module_db, verbose):
    """ Check for any cycles in the complete module dependency graph
    """
    # Build adjacency list
    adjlist = {}

    for m in get_modules(module_db):
        adjlist[m] = list(get_public_deps(m, module_db)
                          | get_private_deps(m, module_db)
                          | get_interface_deps(m, module_db))

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

        print(f"No cycles detected in module dependency graph")
    except graphlib.CycleError as err:
        print(
            f"Cycles detected in module dependency graph:",
            file=sys.stderr,
        )
        cycle = " -> ".join(reversed(err.args[1]))
        print(f"{cycle}", file=sys.stderr)
        sys.exit(1)


def do_find_inconsistencies(db, verbose, full_scan):
    inconsistent = {}
    for m in get_modules(db):
        cc = check_consistency(m, db, full_scan)
        if cc:
            inconsistent[m] = "\n".join(cc)

    if len(inconsistent) > 0:
        for k, v in inconsistent.items():
            print(f'{k}:\n{v}', file=sys.stderr)

        sys.exit(1)
    else:
        print(
            "No inconsistencies detected in declared/apparent module dependencies")


def do_libraries(db, verbose):
    libmap = defaultdict(set)

    for m in get_modules(db):
        libmap[get_parent_target(m, db)].add(m)

    for k, v in sorted(libmap.items()):
        print(f'{k}: {v}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=str(__doc__), formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-db", default='G4ModuleInterfaceMap.csv',
                        metavar="FILE", help="module interface map file")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="verbose output")
    parser.add_argument("--all-files",
                        action="store_true",
                        help="consider all files in source tree for consistency checks")

    query_group = parser.add_mutually_exclusive_group(required=True)

    query_group.add_argument(
        "-l", "--list",
        action='store_true',
        help="list declared source code modules"
    )
    query_group.add_argument(
        "-s", "--source",
        metavar="<module>",
        help="print directory holding CMake file where module is defined"
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
        help="print usage requirements of module"
    )
    query_group.add_argument(
        "-r", "--requires",
        metavar="<module>",
        help="print modules linking to this module"
    )
    query_group.add_argument(
        "--count-usage",
        action="store_true",
        help="print modules from most to least linked to"
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
    query_group.add_argument(
        "--library",
        metavar="<module>",
        help="print final library (.so/.dll) that module is compiled into"
    )
    query_group.add_argument(
        "--libraries",
        action="store_true",
        help="print final libraries and their module compositions"
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
        print("\n".join(get_modules(db)))
    elif args.source:
        try:
            print(get_location(args.source, db))
        except KeyError as err:
            print(f"No module named {err}", file=sys.stderr)
            sys.exit(1)
    elif args.interface:
        try:
            print(
                "\n".join(sorted(get_headers(args.interface, db))))
        except KeyError as err:
            print(f"No module named {err}", file=sys.stderr)
            sys.exit(1)
    elif args.provides:
        do_provides(args.provides, db, args.verbose)
    elif args.usage_requirements:
        ur = usage_requirements(args.usage_requirements, db, args.all_files)
        print(json.dumps(ur, indent=2))
    elif args.count_usage:
        do_count_usage(db, args.verbose)
    elif args.requires:
        do_requires(args.requires, db)
    elif args.check_consistency:
        do_check_consistency(args.check_consistency, db, args.all_files)
    elif args.find_cycles:
        do_find_cycles(db, args.verbose)
    elif args.find_inconsistencies:
        do_find_inconsistencies(db, args.verbose, args.all_files)
    elif args.library:
        try:
            print(get_parent_target(args.library, db))
        except KeyError as err:
            print(f"No module named {err}", file=sys.stderr)
            sys.exit(1)
    elif args.libraries:
        do_libraries(db, args.verbose)
