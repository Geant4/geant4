#! /usr/bin/perl
#
# File: g4rootio.pl
#
#  This script reads the ROOT I/O class defition and implementation
#  information from STDIN or from a *.rootio file.  It then produces the
#  actual C++ header and source files according to the naming rules and
#  directory conventions defined in the style options part of this script.
#
#  To obtain the command list, invoke g4rootio.pl without any option,
#  and type "help".
#
#    % g4rootio.pl
#    fadscpp> help
#    fadscpp> set help
#    fadscpp> quit
#
#  To run in batchmode, specify the name of the file which contains
#  the *.rootio file information.
#
#    % g4rootio.pl [options] <file>
#
#  Option:
#           -h : print usage
#           -f : replace existing source and header files
#           -v : verbose output
#           -n : does not write files
#           -d : difference between new and orginal files
#
# Example:
#  See RPCHit.rootio attached in this distribution for the example
#  *.rootio file.
#
# Style options:
#  By default, above example skeleton will produce the header file
#  as include/MyClassName.hh, source file as src/MyClassName.cc.
#
#  Edit the section "##### Style options #####" below, to change
#  some default behavior.
#
#  Header and source file names can also be specified explicitly
#  with "set header_file <filename>" and "set source_file <filename>".
#
# Caveat:
#  No syntax checking or context checking are done.  For example,
#  puting member declaration into method scope, will create a corrupted
#  source file, but not detected in this perl script.
#  When in doubt, use -v option to see if the skeleton information
#  is correctly loaded.
#
# Author: <youhei.morita@kek.jp>
#

use File::Basename; 
use File::Compare;
use Text::ParseWords;

##### Style options #####

#$author_name       = "enter_your_name";
$author_name        = "";
$package_name       = "Geant4";
$script_name        = basename($0,".pl");
$header_ext_default = ".hh";
$src_ext_default    = ".cc";
$include_dir        = "include";
$src_dir            = "src";
$tmp_dir            = "/tmp";

##### Get the author name from GCOS #####

if ( $author_name eq "" ) {
  chop($user_name = `whoami`);  # chop off the last newline
  ($na,$pw,$uid,$gid,$qo,$cm,$author_name,$dr,$sh) = getpwnam("$user_name");
}

$creation_date=`date +%y.%m.%d`;
chop($creation_date);

##### Parse command line #####

$infile  = "";
$console = 0;
$force   = 0;
$diff    = 0;
$summary = 0;
$nowrite = 0;
$verbose = 0;
$friends[0]      = "";
$global_decls[0] = "";
$assoc_class[0]  = "";
$global_impls[0] = "";
$additional_hdrs_src[0] ="";

if ($#ARGV >= 0 ) {
  while(1) {
    if ( $ARGV[0] =~ /-h/ ) {
      &print_usage;
      exit 1;
    } elsif ( $ARGV[0] =~ /-f/ ) {
      $force = 1;
      shift;
    } elsif ( $ARGV[0] =~ /-d/ ) {
      $diff = 1;
      shift;
    } elsif ( $ARGV[0] =~ /-s/ ) {
      $nowrite = 1;
      $summary = 1;
      shift;
    } elsif ( $ARGV[0] =~ /-n/ ) {
      $nowrite = 1;
      shift;
    } elsif ( $ARGV[0] =~ /-v/ ) {
      $verbose = 1;
      shift;
    } elsif ( $ARGV[0] eq "" ) {
      $infile = "<&STDIN";
      $console = 1;
      last;
    } elsif ( $ARGV[1] eq "" ) {
      $infile = $ARGV[0];
      $force = 1;
      last;
    } else {
      &print_usage;
      exit 1;
    }
  }
} else {
  $infile = "<&STDIN";
  $console = 1;
}

if ( $console == 0 && ! -f $infile ) {
  print "File $infile does not exist.\n";
  exit (1);
}

##### Message #####

if ( $verbose ) {
  print "\n";
  print "***********************************************************\n";
  print "* FADS ROOT I/O adopter class generator                   *\n";
  print "*   utility to generate ROOT I/O header and source files. *\n";
  print "***********************************************************\n\n";
  
  print "File override option is selected.\n" if ( $force );
}

##### Read skeleton information #####

$class_name       = "";
$collection_class = "";
$collection_class_header = 0;
$collection_class_method_entries = "entries()";
$collection_class_method_add     = "insert(obj)";
$collection_class_method_get     = "[i]";
$collection_class_signature = "\@collection_class\@(f_detName, f_colName)";
$collection_base_class = "";
$array_io_base    = "";
$catalog          = "";
$sdet_name        = "";
$make_transient   = "MakeTransientObject";    # method name for retrieve
$template         = "";
$header_file_in   = "";
$source_file_in   = "";
$incl_mark_in     = "";
$namespace        = "none";
$root_class       = 1;
$inherited_class[0] = "TObject";
$num_inherit      = 1;
$root_class_num   = 1;
$num_constructor  = 0;
$num_destructor   = 0;
$num_method       = 0;
$num_member       = 0;
$end              = 0;

print "${script_name}> " if ( ! $force );

open(IN, $infile);

while(<IN>) {
  chop;        # chop off the last newline
  s/^#.*$//;   # omit comments
  @_cmds = quotewords(" ", 1, $_);
  $_cmd = $_cmds[0];
  if ( $_cmd ne "" ) {
    if ( $_cmd =~ /^help$/i ) {
      &fads_help;
    } elsif ( $_cmd =~ /^end$/i ) {
      $end = 1;
      print "'end' accepted.\n" if ( $verbose );
      last;
    } elsif ( $_cmd =~ /^quit$/i ) {
      exit 0;
    } elsif ( $_cmd =~ /^set$/i ) {
      &fads_set;
    } elsif ( $_cmd =~ /^\+EOD$/i ) {
      last;
    } else {
      print "\n-- Unrecognized command: $_cmd\n";
      &fads_help;
      if ( $force ) {
        print "exiting...\n\n";
        exit (1);
      }
    }
  }
  print "${script_name}> " if ( ! $force );
}

print "EOF reached.\n" if ( ! $end && $verbose );

close(IN);

if ( $class_name eq "") {
  print "\n${script_name}: Error - class name is not specified.  Stop.\n";
  exit 1;
}
if ( $collection_class eq "") {
  print "\n${script_name}: Error - collection class is not specified.  Stop.\n";
  exit 1;
}
if ( $collection_base_class eq "") {
  print "\n${script_name}: Error - collection base class is not specified.  Stop.\n";
  exit 1;
}
if ( $array_io_base eq "") {
  print "\n${script_name}: Error - array I/O base class is not specified.  Stop.\n";
  exit 1;
}
if ( $catalog eq "") {
  print "\n${script_name}: Error - catalog is not specified.  Stop.\n";
  exit 1;
}

if ( $sdet_name eq "") {
  print "\n${script_name}: Error - sensitive detctor name is not specified.  Stop.\n";
  exit 1;
}

#----- compose ROOT I/O related class names

$class_root   = $class_name . "Root";
$class_array  = $class_name . "sRoot";
$class_io     = $class_name . "RootIO";
$entry_object = &lowerstring($class_name) . "_" . &lowerstring($catalog);

#----- compose default file names and include header mark
$header_ext = $header_ext_default;

$header_file = "$include_dir/${class_root}${header_ext}";
$source_file = "$src_dir/${class_root}${src_ext_default}";

$header_file = $header_file_in if ( $header_file_in ne "" );
$source_file = $source_file_in if ( $source_file_in ne "" );

#----- temporary files

$header_file_tmp = "$tmp_dir/${class_root}${header_ext}.$$";
$source_file_tmp = "$tmp_dir/${class_root}${src_ext_default}.$$";

#----- ROOT I/O files

$array_header_file = "$include_dir/${class_array}${header_ext}";
$array_source_file = "$src_dir/${class_array}${src_ext_default}";
$io_header_file    = "$include_dir/${class_io}${header_ext}";
$io_source_file    = "$src_dir/${class_io}${src_ext_default}";

$array_header_file_tmp = "$tmp_dir/${class_array}${header_ext}.$$";
$array_source_file_tmp = "$tmp_dir/${class_array}${src_ext_default}.$$";
$io_header_file_tmp    = "$tmp_dir/${class_io}${header_ext}.$$";
$io_source_file_tmp    = "$tmp_dir/${class_io}${src_ext_default}.$$";

#----- Include file mark

$incl_mark = &capitalize(&underscore($class_root));
$incl_mark = $incl_mark. "_HH";
$incl_mark = $incl_mark_in if ( $incl_mark_in ne "" );

#----- Convert macros

&fads_conv;

#----- Print information

if ( $verbose ) {
  print "Include Header Mark: $incl_mark\n";
  print "Default Constructor will be used.\n" if ( $constructors[0] eq "" );
  print "Default Destructor will be used.\n" if ( $destructors[0] eq "" );
}

#----- Check existence of header and source directories

if ( ! $nowrite ) {
  &check_dir($include_dir);
  &check_dir($src_dir);
}

##### Create tmp header file #####

&fads_create_header_file;

##### Create tmp source file #####

$no_source = &is_no_source;
if ( ! $no_source ) {
  &fads_create_source_file;
} else {
  print "Warning: file $source_file is _not_ created.\n";
  print " - All methods are inline, pure virtual, or the class is template.\n";
}

##### Create tmp ROOT I/O files #####

&fads_create_array_header_file;
&fads_create_array_source_file;
&fads_create_io_header_file;
&fads_create_io_source_file;

#----- diff-ing files for -n option

if ( $diff ) {
  $_ = 1;                               # does fflush(3) for each print

  &diff_file( $header_file_tmp, $header_file );
  &diff_file( $source_file_tmp, $source_file ) if ( ! $no_source );
  &diff_file( $array_header_file_tmp, $array_header_file );
  &diff_file( $array_source_file_tmp, $array_source_file );
  &diff_file( $io_header_file_tmp, $io_header_file );
  &diff_file( $io_source_file_tmp, $io_source_file );
}

#----- summary output for -s option

if ( $summary ) {
  &fads_print_summary;
}

#----- Install header and source files

if ( ! $nowrite ) {
  &check_file($header_file);
  &check_file($source_file) if ( ! $no_source );
  &check_file($array_header_file);
  &check_file($array_source_file);
  &check_file($io_header_file);
  &check_file($io_source_file);

  &install_file( $header_file_tmp, $header_file );
  if ( ! $no_source ) {
    &install_file( $source_file_tmp, $source_file );
  } else {
    system("rm -f $source_file") if ( -f $source_file );
  }
  &install_file( $array_header_file_tmp, $array_header_file );
  &install_file( $array_source_file_tmp, $array_source_file );
  &install_file( $io_header_file_tmp, $io_header_file );
  &install_file( $io_source_file_tmp, $io_source_file );
}

#----- Delete tmp header and source files

system("rm -f $header_file_tmp");
system("rm -f $source_file_tmp") if ( ! $no_source );
system("rm -f $array_header_file_tmp");
system("rm -f $array_source_file_tmp");
system("rm -f $io_header_file_tmp");
system("rm -f $io_source_file_tmp");

#----- End of the script

exit 0;

############################################################################

sub print_usage {
  print "\n Usage:  ${script_name} [-h] [-f] [-d] [-s] [-n] [-v] <file>\n\n";
  print "     -h : print usage\n";
  print "     -f : override existing header and source files\n";
  print "     -d : difference output to stdout\n";
  print "     -s : summary output to stdout\n";
  print "     -n : does not replace files\n";
  print "     -v : verbose mode\n";
  print "\n";
}

############################################################################

sub fads_help {
  print "\n Available commands:\n";
  print "    help   :  prints this message\n";
  print "    set    :  sets various parameters\n";
  print "    end    :  ends class description and creates skeleton files\n";
  print "    quit   :  quits this program without changing files\n";
  print "\n";
}

############################################################################

sub fads_set_help {
  print "\n Available set command options:\n";
  print "  help                  : prints this message\n";
  print "  class_name            : sets the class name\n";
  print "  collection_class      : sets the collection class name\n";
  print "  collection_class_header : set yes if the collection has a header file\n";
  print "  collection_class_method_entries :   method name for entries\n";
  print "  collection_class_method_add     :   method name to add an entry\n";
  print "  collection_class_method_get     :   method name to get i-th entry\n";
  print "  collection_class_signagure      :   constructor signature\n";
  print "  collection_base_class : sets the collection base class name\n";
  print "  array_io_base         : sets the array I/O base class name\n";
  print "  catalog               : sets the catalog name\n";
  print "  sdet_name             : sets the sensitive detector name\n";
  print "  template              : sets the class template\n";
  print "  description           : sets the class description\n";
  print "  history               : sets the development history\n";
  print "  header_file           : sets the header file name\n";
  print "  source_file           : sets the source file name\n";
  print "  include_mark          : sets the include macro name\n";
  print "  namespace             : sets the name space\n";
  print "  root_class_num        : sets the root class number\n";
  print "  inherited_class       : class to be inherited\n";
  print "  associated_class      : associated class to be used\n";
  print "  add_header            : additional header include\n";
  print "  add_header_src        : additional header include for source file\n";
  print "  global_declaration    : add. global declaration in header file\n";
  print "  global_implementation : add. global impl. in source file\n";
  print "  friend                : sets the friend class declaration\n";
  print "  constructor           : explicitly create constructor(s)\n";
  print "  destructor            : explicitly create destructor(s)\n";
  print "  method [ public|protected|private ] \n";
  print "                        : set method types, names and arguments\n";
  print "  member [ public|protectec|private ] \n";
  print "                        : set data member types and names\n";
  print "\n";
}

############################################################################

sub fads_set {
  $option = $_cmds[1];
  if ( $option =~ /^help$/i ) {
          &fads_set_help;
  } elsif ( $option =~ /^class_name$/i ) {
               &fads_set_class_name;
  } elsif ( $option =~ /^collection_class$/i ) {
               &fads_set_collection_class;
  } elsif ( $option =~ /^collection_class_header$/i ) {
               &fads_set_collection_class_header;
  } elsif ( $option =~ /^collection_class_method_entries$/i ) {
               &fads_set_collection_class_method_entries;
  } elsif ( $option =~ /^collection_class_method_add$/i ) {
               &fads_set_collection_class_method_add;
  } elsif ( $option =~ /^collection_class_method_get$/i ) {
               &fads_set_collection_class_method_get;
  } elsif ( $option =~ /^collection_class_signature$/i ) {
               &fads_set_collection_class_signature;
  } elsif ( $option =~ /^collection_base_class$/i ) {
               &fads_set_collection_base_class;
  } elsif ( $option =~ /^array_io_base$/i ) {
               &fads_set_array_io_base;
  } elsif ( $option =~ /^catalog$/i ) {
               &fads_set_catalog;
  } elsif ( $option =~ /^sdet_name$/i ) {
               &fads_set_sdet_name;
  } elsif ( $option =~ /^template$/i ) {
               &fads_set_template;
  } elsif ( $option =~ /^description$/i ) {
               &fads_set_description;
  } elsif ( $option =~ /^history$/i ) {
               &fads_set_history;
  } elsif ( $option =~ /^header_file$/i ) {
               &fads_set_header_file;
  } elsif ( $option =~ /^source_file$/i ) {
               &fads_set_source_file;
  } elsif ( $option =~ /^include_mark$/i ) {
               &fads_set_include_mark;
  } elsif ( $option =~ /^namespace$/i ) {
               &fads_set_namespace;
  } elsif ( $option =~ /^root_class_num$/i ) {
               &fads_set_root_class_num;
  } elsif ( $option =~ /^friend$/i ) {
               &fads_set_friend;
  } elsif ( $option =~ /^constructor$/i ) {
               &fads_set_constructor;
  } elsif ( $option =~ /^destructor$/i ) {
               &fads_set_destructor;
  } elsif ( $option =~ /^inherited_class$/i ) {
               &fads_set_inherited_class;
  } elsif ( $option =~ /^associated_class$/i ) {
               &fads_set_associated_class;
  } elsif ( $option =~ /^add_header$/i ) {
               &fads_set_add_header;
  } elsif ( $option =~ /^add_header_src$/i ) {
               &fads_set_add_header_src;
  } elsif ( $option =~ /^global_declaration$/i ) {
               &fads_set_global_declaration;
  } elsif ( $option =~ /^global_implementation$/i ) {
               &fads_set_global_implementation;
  } elsif ( $option =~ /^method$/i ) {
               &fads_set_method;
  } elsif ( $option =~ /^member$/i ) {
               &fads_set_members;
  } else {
    print "\n${script_name}: Unrecognized set command option: $option\n";
    &fads_set_help;
    if ( $force ) {
      print "exiting...\n\n";
      exit (1);
    }
  }
}

############################################################################

sub fads_set_class_name {
  if ( $#_cmds < 2 ) {
    print "\nError - specify a class name.\n";
    exit(1) if ( $force );
    return;
  }
  $class_name = "$_cmds[2]";
  print "Class name = $class_name\n" if ($verbose) ;
}

############################################################################

sub fads_set_collection_class {
  if ( $#_cmds < 2 ) {
    print "\nError - specify a collection class name.\n";
    exit(1) if ( $force );
    return;
  }
  $collection_class = "$_cmds[2]";
  print "Collection class name = $collection_class\n" if ($verbose) ;
}

############################################################################

sub fads_set_collection_class_header {
  if ( $#_cmds < 2 ) {
    print "\nError - specify yes or no!.\n";
    exit(1) if ( $force );
    return;
  }
  &tidy_command_option;
  if ( $_cmds[2] =~ /^y$/i || $_cmds[2] =~ /^yes$/i  ) {
    $collection_class_header = 1
  } elsif ( $_cmds[2] =~ /^n$/i || $_cmds[2] =~ /^no$/i  ) {
    $collection_class_header = 0
  } else {
    print "\nError - specify yes or no!.\n";
  }
}

############################################################################

sub fads_set_collection_class_method_entries {
  if ( $#_cmds < 2 ) {
    print "\nError - specify a collection class method name for entries.\n";
    exit(1) if ( $force );
    return;
  }
  &tidy_command_option;
  $collection_class_method_entries = "$_cmds[2]";
  print "Collection class method for entries = $collection_class_method_entries\n" if ($verbose) ;
}

############################################################################

sub fads_set_collection_class_method_add {
  if ( $#_cmds < 2 ) {
    print "\nError - specify a collection class method name for add.\n";
    exit(1) if ( $force );
    return;
  }
  &tidy_command_option;
  $collection_class_method_add = "$_cmds[2]";
  print "Collection class method for add = $collection_class_method_add\n" if ($verbose) ;
}

############################################################################

sub fads_set_collection_class_method_get {
  if ( $#_cmds < 2 ) {
    print "\nError - specify a collection class method name for get.\n";
    exit(1) if ( $force );
    return;
  }
  &tidy_command_option;
  $collection_class_method_get = "$_cmds[2]";
  print "Collection class method for get = $collection_class_method_get\n" if ($verbose) ;
}

############################################################################

sub fads_set_collection_class_signature {
  if ( $#_cmds < 2 ) {
    print "\nError - specify a collection class signature.\n";
    exit(1) if ( $force );
    return;
  }
  &tidy_command_option;
  $collection_class_signature = "$_cmds[2]";
  print "Collection class signature = $collection_class_signature\n" if ($verbose) ;
}

############################################################################

sub fads_set_collection_base_class {
  if ( $#_cmds < 2 ) {
    print "\nError - specify a collection base class name.\n";
    exit(1) if ( $force );
    return;
  }
  $collection_base_class = "$_cmds[2]";
  print "Collection base class name = $collection_base_class\n" if ($verbose) ;
}

############################################################################

sub fads_set_array_io_base {
  if ( $#_cmds < 2 ) {
    print "\nError - specify an array I/O base class name.\n";
    exit(1) if ( $force );
    return;
  }
  $array_io_base = "$_cmds[2]";
  print "Array I/O base class name = $array_io_base\n" if ($verbose) ;
}

############################################################################

sub fads_set_catalog {
  if ( $#_cmds < 2 ) {
    print "\nError - specify a catalog name.\n";
    exit(1) if ( $force );
    return;
  }
  $catalog = "$_cmds[2]";
  print "Catalog name = $catalog\n" if ($verbose) ;
}

############################################################################

sub fads_set_sdet_name {
  if ( $#_cmds < 2 ) {
    print "\nError - specify a sensitive detector name.\n";
    exit(1) if ( $force );
    return;
  }
  $sdet_name = "$_cmds[2]";
  print "Sensitive detector name = $sdet_name\n" if ($verbose) ;
}

############################################################################

sub fads_set_template {
  if ( $#_cmds < 2 ) {
    print "\nError - specify template.\n";
    exit(1) if ( $force );
    return;
  }
  shift @_cmds;
  shift @_cmds;
  while ($_cmds[0] ne "") {
    $template .= " " . $_cmds[0];
    shift @_cmds;
  }
  $template =~ s/^ //;
  print "Template = $template\n" if ($verbose) ;
}

############################################################################

sub fads_set_description {
  $i=0;
  while(<IN>) {
    chop;         # chop off the last newline
    s/^#.*$//;    # omit comments
    last if ( &is_endline );
    $description[$i++] = $_;
  }
  if ($verbose) {
    print "--- Class Description: ------------------------------------\n";
    foreach $line (@description) {
      print $line,"\n";
    }
    print "--- Class Description: ------------------------------------\n\n";
  }
}

############################################################################

sub fads_set_history {
  $i=0;
  while(<IN>) {
    chop;         # chop off the last newline
    s/^#.*$//;    # omit comments
    last if ( &is_endline );
    $history[$i++] = $_;
  }
  if ($verbose) {
    print "--- History: ----------------------------------------------\n";
    foreach $line (@history) {
      print $line,"\n";
    }
    print "--- History: ----------------------------------------------\n\n";
  }
}

############################################################################

sub fads_set_header_file {
  if ( $#_cmds < 2 ) {
    print "\nError - specify header file name.\n";
    exit(1) if ( $force );
    return;
  }
  $header_file_in = "$_cmds[2]";
  print "Header file = $header_file_in\n" if ($verbose);
}

############################################################################

sub fads_set_source_file {
  if ( $#_cmds < 2 ) {
    print "\nError - specify source file name.\n";
    exit(1) if ( $force );
    return;
  }
  $source_file_in = "$_cmds[2]";
  print "Source file = $source_file_in\n" if ($verbose);
}

############################################################################

sub fads_set_include_mark {
  if ( $#_cmds < 2 ) {
    print "\nError - specify include mark name.\n";
    exit(1) if ( $force );
    return;
  }
  $incl_mark_in = "$_cmds[2]";
  print "Include mark = $incl_mark_in\n" if ($verbose);
}

############################################################################

sub fads_set_namespace {
  if ( $#_cmds < 2 ) {
    print "\nError - specify name space.\n";
    exit(1) if ( $force );
    return;
  }
  $namespace = "$_cmds[2]";
  print "Name space = $namespace\n" if ($verbose);
}

############################################################################

sub fads_set_root_class_num {
  if ( $#_cmds < 2 ) {
    print "\nError - specify the class number!.\n";
    exit(1) if ( $force );
    return;
  }
  $root_class_num = $_cmds[2];
  if ($verbose) {
    print "*** ROOT class number $root_class_num ***\n"
  }
}

############################################################################

sub fads_set_constructor {
  $num_constructor++;
  $i = $num_constructor - 1;
  $num_constructor_impl[$i] = 0;
  $constructor_desc[$i] = "";
  $constructor_impl[$i] = "";
  $l = 0;
  while(<IN>) {
    chop($line0 = $_);
    $_ = &tidy_line;
    last if ( &is_endline );
    if ( $_ ne "" ) {
      $l++;
      if ( $l == 1 ) {
        ($idt = $line0) =~ s/(^[ ]*)(.*$)/$1/;  # count indent space
        s/;[ ]*$//;
        $constructors[$i] = $_;
      } elsif ( $num_constructor_impl[$i] == 0 && /^\/\// ) {
        s/^\/\/ //;                             # chop off C++ comment mark
        $constructor_desc[$i] .= $_ . "\n";
      } else {
        $line0 =~ s/(^$idt)//;                  # chop off indent space
        $num_constructor_impl[$i]++;
        $constructor_impl[$i] .= $line0 . "\n";
      }
    } else {
      if ( $num_constructor_impl[$i] > 0 ) {
        $num_constructor_impl[$i]++;
        $constructor_impl[$i] .= "\n";
      }
    }
  }
  &fads_print_constructor($i, $constructors[$i]) if ($verbose);
}

############################################################################

sub fads_print_summary {
  print "${class_root}:\n";
  if ( $inherited_class[0] ne "" ) {
    print "  -> ";
    foreach $iclass (@inherited_class) {
      print "${iclass} ";
    }
    print "\n";
  }
  $i=0;
  foreach $line (@member_names) {
    &fads_print_member($i++, "summary");
  }
  $i=0;
  foreach $line (@method_names) {
    &fads_print_method($i++, $line, "summary");
  }
  print " collection_class: $collection_class\n";
  print " collection_class_header: $collection_class_header\n";
  print " collection_class_method_entries: $collection_class_method_entries\n";
  print " collection_class_method_add: $collection_class_method_add\n";
  print " collection_class_method_get: $collection_class_method_get\n";
  print " collection_class_signature:  $collection_class_signature\n";
  print " collection_base_class: $collection_base_class\n";
  print " array_io_base: $array_io_base\n";
  print " catalog: $catalog\n";
}

############################################################################

sub fads_print_constructors {
  $i=0;
  foreach $line (@constructors) {
    &fads_print_constructor($i++, $line);
  }
}

############################################################################

sub fads_print_constructor {
  $j=$_[0]+1;
  print "Constructor #${j}: $_[1]\n";
  &fads_print_description($constructor_desc[$_[0]]);
  &fads_print_implementation_body($constructor_impl[$_[0]],
                                  $num_constructor_impl[$_[0]]);
}

############################################################################

sub fads_set_destructor {
  $num_destructor++;
  $i = $num_destructor - 1;
  $num_destructor_impl[$i] = 0;
  $destructor_desc[$i] = "";
  $destructor_impl[$i] = "";
  $l = 0;
  while(<IN>) {
    chop($line0 = $_);
    $_ = &tidy_line;
    last if ( &is_endline );
    if ( $_ ne "" ) {
      $l++;
      if ( $l == 1 ) {
        ($idt = $line0) =~ s/(^[ ]*)(.*$)/$1/;  # count indent space
        s/;[ ]*$//;
        $destructors[$i] = $_;
      } elsif ( $num_destructor_impl[$i] == 0 && /^\/\// ) {
        s/^\/\/ //;                             # chop off C++ comment mark
        $destructor_desc[$i] .= $_ . "\n";
      } else {
        $line0 =~ s/(^$idt)//;                  # chop off indent space
        $num_destructor_impl[$i]++;
        $destructor_impl[$i] .= $line0 . "\n";
      }
    } else {
      if ( $num_destructor_impl[$i] > 0 ) {
        $num_destructor_impl[$i]++;
        $destructor_impl[$i] .= "\n";
      }
    }
  }
  &fads_print_destructor($i, $destructors[$i]) if ($verbose);
}

############################################################################

sub fads_print_destructors {
  $i=0;
  foreach $line (@destructors) {
    &fads_print_destructor($i++, $line);
  }
}

############################################################################

sub fads_print_destructor {
  $j=$_[0]+1;
  print "Destructor #${j}: $_[1]\n";
  &fads_print_description($destructor_desc[$_[0]]);
  &fads_print_implementation_body($destructor_impl[$_[0]],
                                  $num_destructor_impl[$_[0]]);
}

############################################################################

sub fads_set_friend {
  $i=0;
  while(<IN>) {
    $_ = &tidy_line;
    last if ( &is_endline );
    $friends[$i++] = $_;
  }
  if ($verbose) {
    print "Friend classes: @friends\n";
  }
}

############################################################################

sub fads_set_inherited_class {
  $i=$num_inherit;
  while(<IN>) {
    $_ = &tidy_line;
    last if ( &is_endline );
    $root_class = 1 if ($_ eq "TObject");
    $inherited_class[$i++] = $_;
  }
  $num_inherit = $i;
  if ($verbose) {
    print "Number of inherited class: $num_inherit\n";
    print "Inherited class: @inherited_class\n";
    print "*** ROOT class ***\n" if ($root_class == 1)
  }
}

############################################################################

sub fads_set_associated_class {
  $i=0;
  while(<IN>) {
    $_ = &tidy_line;
    last if ( &is_endline );
    $assoc_class[$i++] = $_;
  }
  print "Associated class: @assoc_class\n" if ($verbose);
}

############################################################################

sub fads_set_add_header {
  $i=0;
  while(<IN>) {
    $_ = &tidy_line;
    last if ( &is_endline );
    $additional_hdrs[$i++] = $_;
  }
  print "Additional header include: @additional_hdrs\n" if ($verbose);
}

############################################################################

sub fads_set_add_header_src {
  $i=0;
  while(<IN>) {
    $_ = &tidy_line;
    last if ( &is_endline );
    $additional_hdrs_src[$i++] = $_;
  }
  print "Additional header include for source file: @additional_hdrs_src\n"
                                                              if ($verbose);
}

############################################################################

sub fads_set_global_declaration {
  $i=0;
  while(<IN>) {
    $_ = &tidy_line;
    last if ( &is_endline );
    $global_decls[$i++] = $_;
  }
  print "Global declation: @global_decls\n" if ($verbose);
}

############################################################################

sub fads_set_global_implementation {
  $i=0;
  while(<IN>) {
    $_ = &tidy_line;
    last if ( &is_endline );
    $global_impls[$i++] = $_;
  }
  print "Global implementation: @global_impls\n" if ($verbose);
}

############################################################################

sub fads_set_method {
  $num_method++;
  $i = $num_method - 1;
  $num_method_impl[$i] = 0;
  $method_desc[$i] = "";
  $method_impl[$i] = "";
  $method_attr[$i] = "public";
  $method_attr[$i] = "protected" if ($_cmds[2] eq "protected");
  $method_attr[$i] = "private"   if ($_cmds[2] eq "private");
  $l = 0;
  while(<IN>) {
    chop($line0 = $_);
    $_ = &tidy_line;
    last if ( &is_endline );
    if ( $_ ne "" ) {
      $l++;
      if ( $l == 1 ) {
        ($idt = $line0) =~ s/(^[ ]*)(.*$)/$1/;  # count indent space
        s/;[ ]*$//;                             # chop off trailing semicolon
        $method_names[$i] = $_;
      } elsif ( $num_method_impl[$i] == 0 && /^\/\// ) {
        s/^\/\/ //;                             # chop off C++ comment mark
        $method_desc[$i] .= $_ . "\n";
      } else {
        $line0 =~ s/(^$idt)//;                  # chop off indent space
        if ( &test_line($line0) ) {
          $num_method_impl[$i]++;
          $method_impl[$i] .= $line0 . "\n";
        } else {
          print "\nError: Line does not match.  Check the number of double quotation.  Stop.\n";
          print "  Class:  $class_name\n";
          print "  Method: $method_names[$i]\n";
          print "  Line:   $line0\n";
          exit 1;
        }
      }
    } else {
      if ( $num_method_impl[$i] > 0 ) {
        $num_method_impl[$i]++;
        $method_impl[$i] .= "\n";
      }
    }
  }
  &fads_print_method($i, $method_names[$i]) if ( $verbose );
}

############################################################################

sub test_line {
  $_ = $_[0] . "\n";
  @words = quotewords("\n",1,$_);
  if ( $words[0] ne "" ) {
    return 1;
  } else {
    return 0;
  }
}

############################################################################

sub fads_print_methods {
  $i=0;
  foreach $line (@method_names) {
    &fads_print_method($i++, $line);
  }
}

############################################################################

sub fads_print_method {
  $j=$_[0]+1;
  if ( $_[2] eq "summary" ) {
    $rtype  = &return_type($_[1]);
    $namarg = &method_namarg($_[1]);
    $namarg =~ s/(^.*)\{.*$/$1/;
    print "  + $namarg  :  $rtype\n";
  } else {
    print "Method#${j} ";
    print "$method_attr[$_[0]]: ";
    print "virtual " if ( &is_virtual($_[1]) );
    print "static "  if ( &is_static($_[1]) );
    print "inline "  if ( &is_inline($_[1]) );
    $rtype  = &return_type($_[1]);
    $namarg = &method_namarg($_[1]);
    print "$rtype $namarg\n";
    &fads_print_description($method_desc[$_[0]]);
    &fads_print_implementation_body($method_impl[$_[0]],
                                    $num_method_impl[$_[0]] );
  }
}

############################################################################

sub fads_set_members {
  $attr = "private";
  $attr = "protected" if ($_cmds[2] eq "protected");
  $attr = "public"    if ($_cmds[2] eq "public");
  while(<IN>) {
    $_ = &tidy_line();
    last if ( &is_endline() );
    if ( $_ ne "" ) {
      &fads_set_member($num_member++, $_, $attr);
    }
  }
}

############################################################################

sub fads_set_member {
  $_ = $_[1];
  s/;[ ]*$//;
  $type = $_;
  $type =~ s/(^.*)[ ]+(.*$)/$1/g;
  $mname = $_;
  $mname =~ s/(^.*)[ ]+(.*$)/$2/g;
  $member_type[$_[0]]  = $type;
  $member_names[$_[0]] = $mname;
  $member_attr[$_[0]]  = $_[2];
  $member_attr[$_[0]]  = "private" if ( $_[2] eq "" );

  &fads_print_member($_[0]) if ( $verbose );
}

############################################################################

sub fads_print_members {
  $i=0;
  foreach $line (@member_names) {
    &fads_print_member($i++);
  }
}

############################################################################

sub fads_print_member {
  $j=$_[0]+1;
  if ( $_[1] eq "summary" ) {
    print "  - $member_type[$_[0]] $member_names[$_[0]]\n";
  } else {
    print "member#${j} $member_attr[$_[0]]: ";
    print "$member_type[$_[0]] $member_names[$_[0]]\n";
  }
}

############################################################################

sub fads_conv {

  $collection_class_method_entries = &macro_conv($collection_class_method_entries);

  # change method_add arguments to "(a)"
  $collection_class_method_add =~ s/(.*)\(.*\)/$1(a)/;
  $collection_class_method_add = &macro_conv($collection_class_method_add);

  # add "." to method_get.  Don't add if [] operator.
  $collection_class_method_get =~ s/^(.*)$/.$1/;
  $collection_class_method_get =~ s/^.\[(.*)$/[$1/;
  $collection_class_method_get = &macro_conv($collection_class_method_get);

  $collection_class_signature = &macro_conv($collection_class_signature);

  $i=0;
  foreach $line (@additional_hdrs) {
    $additional_hdrs[$i] = &macro_conv($additional_hdrs[$i]);
    $i++;
  }

  $i=0;
  foreach $line (@additional_hdrs_src) {
    $additional_hdrs_src[$i] = &macro_conv($additional_hdrs_src[$i]);
    $i++;
  }

  $i=0;
  foreach $line (@global_decls) {
    $global_decls[$i] = &macro_conv($global_decls[$i]);
    $i++;
  }

  $i=0;
  foreach $line (@global_impls) {
    $global_impls[$i] = &macro_conv($global_impls[$i]);
    $i++;
  }

  $i=0;
  foreach $line (@member_names) {
    $member_type [$i] = &macro_conv($member_type [$i]);
    $member_names[$i] = &macro_conv($member_names[$i]);
    $i++;
  }

  $i=0;
  foreach $line (@method_names) {
    $method_desc [$i] = &macro_conv($method_desc [$i]);
    $method_impl [$i] = &macro_conv($method_impl [$i]);
    $method_names[$i] = &macro_conv($method_names[$i]);
  }

  $i=0;
  foreach $line (@constructors) {
    $constructor_desc[$i] = &macro_conv($constructor_desc[$i]);
    $constructor_impl[$i] = &macro_conv($constructor_impl[$i]);
    $constructors[$i]     = &macro_conv($constructors[$i]);
  }

  $i=0;
  foreach $line (@destructor_names) {
    $destructor_desc[$i] = &macro_conv($destructor_desc [$i]);
    $destructor_impl[$i] = &macro_conv($destructor_impl [$i]);
    $destructors[$i]     = &macro_conv($destructors[$i]);
  }
}

############################################################################

sub fads_print_implementation_body {
  @words = quotewords("\n",1,$_[0]);
  $k=0;
  foreach $word (@words) {
    $k++;
    print "  $word\n" if ( $k <= $_[1] );
  }
}

############################################################################

sub macro_conv {
  $_ = $_[0];
  s/\@class_name\@/$class_name/g;
  s/\@class_root\@/$class_root/g;
  s/\@class_array\@/$class_array/g;
  s/\@class_io\@/$class_io/g;
  s/\@collection_class\@/$collection_class/g;
  s/\@collection_class_method_entries\@/$collection_class_method_entries/g;
  s/\@collection_class_method_add\@/$collection_class_method_add/g;
  s/\@collection_class_method_get\@/$collection_class_method_get/g;
  s/\@collection_class_signature\@/$collection_class_signature/g;
  s/\@collection_base_class\@/$collection_base_class/g;
  s/\@array_io_base\@/$array_io_base/g;
  s/\@catalog\@/$catalog/g;
  s/\@entry_object\@/$entry_object/g;
  s/\@sdet_name\@/$sdet_name/g;
  s/\@make_transient\@/$make_transient/g;

  # for ROOT data types
  s/\@int\@/ Int_t/g;
  s/\@uint\@/UInt_t/g;
  s/\@float\@/ Float_t/g;
  s/\@char\@/char/g;

  return $_;
}

############################################################################

sub capitalize {
  $_ = $_[0];
  y/a-z/A-Z/;
  return $_;
}

############################################################################

sub lowerstring {
  $_ = $_[0];
  y/A-Z/a-z/;
  return $_;
}

############################################################################

sub underscore {
  $_ = $_[0];
  s/[A-Z]/_$&/g;
  s/^_//g;
  s/_([A-Z])_([A-Z])_/$1$2_/g;
  s/_([A-Z])_([A-Z])$/_$1$2/g;
  return $_;
}

############################################################################

sub tidy_command_option {
  $i = 0;
  foreach $cmd (@_cmds) {
    last if ( $cmd eq "#" );
    if ( $i > 2 ) {
      $_tmp[2] .= $cmd if ( $cmd ne "" );
    } else {
      $_tmp[$i++] = $cmd if ( $cmd ne "" );
    }
  }
  @_cmds = ();
  while( $_ = shift(@_tmp) ) {
    push(@_cmds,$_);
  }
}

############################################################################

sub tidy_line {
  $result = $_[0];
  $result = $_ if ( $result eq "" );
  chop($result);         # chop off the last newline
  $result =~ s/^#.*$//;  # omit comments
  $result =~ s/^[ ]+//;  # chop off leading spaces
  $result =~ s/[ ]+$//;  # chop off the trailing branks
##  $result =~ s/;.*$//;   # chop off the rest of the words after ";"
  return $result;
}

############################################################################

sub fads_create_header_file {
  open(OUT, "> $header_file_tmp");

  print(OUT "//\n");
  print(OUT "// File: ${class_root}${header_ext}\n");
  print(OUT "//\n");
  print(OUT "//   This file is automatically generated by $script_name.pl.\n");
  print(OUT "//   Do not edit by hand.  Edit $class_name.rootio instead.\n");
  print(OUT "\n");

  print(OUT "#ifndef $incl_mark\n");
  print(OUT "#define $incl_mark 1\n");
  print(OUT "\n");

  if ( $additional_hdrs[0] ne "" ) {
    foreach $addhdr (@additional_hdrs) {
      $line = "#include \"$addhdr\"";
      $line =~ s/"</</;
      $line =~ s/>"/>/;
      print(OUT "$line\n");
    }
    print(OUT "\n");
  }

  if ( $global_decls[0] ne "" ) {
    foreach $decl (@global_decls) {
      print(OUT "$decl\n");
    }
    print(OUT "\n");
  }

  if ( $inherited_class[0] ne "" ) {
##    print(OUT "// Class inherited:\n");
    foreach $iclass (@inherited_class) {
      if ( $iclass ne "TObject" ) {
        print(OUT "#include \"${iclass}${header_ext_default}\"\n");
      } else {
        print(OUT "#include \"TObject.h\"\n");
      }
    }
    print(OUT "\n");
  }

  if ( $assoc_class[0] ne "" ) {
    print(OUT "// Forward Declarations:\n");
    foreach $aclass (@assoc_class) {
      print(OUT "class ${aclass};\n");
    }
    print(OUT "\n");
  }

  if ( $description[0] ne "" ) {
    print(OUT "// Class Description:\n");
    foreach $line (@description) {
      print(OUT "// $line\n");
    }
    print(OUT "\n");
  }

  if ( $namespace ne "none" ) {
    print(OUT "namespace ${namespace}\n");
    print(OUT "{\n");
    print(OUT "\n");
    $nspace = "  ";
  } else {
    $nspace = "";
  }

  if ( $template ne "" ) {
    print(OUT "template $template ");
  }

  print(OUT "${nspace}class $class_root\n");

  if ( $inherited_class[0] ne "" ) {
    print(OUT "${nspace} : ");
    $i=0;
    foreach $iclass (@inherited_class) {
      $i++;
      print(OUT "public $iclass");
      print(OUT ", ")  if ( $i != $num_inherit );
    }
  }
  print(OUT "\n") if ( $num_inherit > 0 );

  print(OUT "${nspace}\{\n");

  if ( $friends[0] ne "" ) {
    foreach $friend (@friends) {
      print(OUT "    friend class ${friend};\n") if ( $friend ne "" );
    }
    print(OUT "\n");
  }

  print(OUT "    public:\n");
  if ( $constructors[0] ne "" ) {
    $i = 0;
    foreach $constructor (@constructors) {
      if ( $constructor ne "" ) {
        print(OUT "      ${constructor}");
        if ( $template ne "" ) {
          if ( $num_constructor_impl[$i] > 0 ) {
            print(OUT "\n")
            &fads_create_implementation_body($constructor_impl[$i],
                                         $num_constructor_impl[$i], "      ");
          } elsif ( ! &is_implicit_inline($constructor) ) {
            print(OUT "{};\n");
          }
        } else {
          print(OUT ";\n")
        }
      }
      $i++;
    }
  } else {
    if ( $template ne "" ) {
      print(OUT "      ${class_root}() {};\n");
    } else {
      print(OUT "      ${class_root}();\n");
    }
  }
##  &fads_create_description($constructor_desc[$i], "Constructor");

  if ( $destructors[0] ne "" ) {
    foreach $destructor (@destructors) {
      print(OUT "      ${destructor};\n") if ( $destructor ne "" );
    }
  } else {
    print(OUT "      ~${class_root}();\n");
  }
##  &fads_create_description($destructor_desc[$i], "Destructor");
  print(OUT "\n");

  &fads_create_header_methods("public");
  &fads_create_header_methods("protected");
  &fads_create_header_methods("private");

  if ( $assoc_class[0] ne "" ) {
    print(OUT "    private:\n");
    foreach $aclass (@assoc_class) {
      print(OUT "      ${aclass}* f_${aclass};\n");
    }
    print(OUT "\n");
  }

  &fads_create_header_members("public");
  &fads_create_header_members("protected");
  &fads_create_header_members("private");

  if ( $root_class == 1 ) {
    print(OUT "    ClassDef (${class_root},${root_class_num})\n");
    print(OUT "\n");
  }

  print(OUT "${nspace}}; // End of class ${class_root}\n");
  print(OUT "\n");

  if ( $namespace ne "none" ) {
    print(OUT "} // ${namespace}\n");
    print(OUT "\n");
  }

  print(OUT "#endif\n");
  print(OUT "\n");

  close(OUT);

  print "File $header_file created.\n" if ( $verbose && ! $nowrite );
}

############################################################################

sub fads_create_header_methods {
  $attr = $_[0];
  $attr = "public" if ( $_[0] eq "");
  if ( &check_methods_attr($attr) ) {
    print(OUT "    $attr:\n");
    $i=0;
    foreach $method (@method_names) {
      if ( $method ne "" ) {
        &fads_create_header_method($i++, $method, $attr)
      }
    }
    print(OUT "\n");
  }
}

############################################################################

sub fads_create_header_method {
  if ( $_[2] eq $method_attr[$_[0]] ) {
    $j = $_[0]+1;
    print(OUT "      ");
    print(OUT "inline ")  if ( &is_inline($_[1]) );
    print(OUT "static ")  if ( &is_static($_[1]) );
    print(OUT "virtual ") if ( &is_virtual($_[1]) );
    $rtype  = &return_type($_[1]);
    $rtype =~ s/\<namespace\>//;
    $namarg = &method_namarg($_[1]);
    print(OUT "$rtype $namarg");
    if ( $template ne "" ) {
      if ( $num_method_impl[$_[0]] > 0 ) {
        print(OUT "\n");
        @words = quotewords("\n",1,$method_impl[$_[0]]);
        $k=0;
        foreach $word (@words) {
          $k++;
          print(OUT "      $word\n") if ( $k <= $num_method_impl[$_[0]] );
        }
      } else {
        print(OUT ";\n");
      }
    } else {
      print(OUT ";\n");
    }
    $name = "method #$j: $namarg";
##    &fads_create_description($method_desc[$_[0]], $name);
  }
}

############################################################################

sub fads_create_header_members {
  $attr = $_[0];
  $attr = "private" if ( $_[0] eq "");
  if ( &check_members_attr($attr) ) {
    print(OUT "    $attr:\n");
    $i=0;
    foreach $member (@member_names) {
      if ( $member ne "" ) {
        &fads_create_header_member($i++, $attr);
      }
    }
    print(OUT "\n");
  }
}

############################################################################

sub fads_create_header_member {
  $mtype = $member_type[$_[0]]; 
  $mtype =~ s/\<namespace\>/${namespace}::/  if ( $namespace ne "none" );
  $mtype =~ s/\<namespace\>//                if ( $namespace eq "none" );
  if ( $_[1] eq $member_attr[$_[0]] ) {
    print(OUT "      $mtype $member_names[$_[0]];\n");
  }
}

############################################################################

sub check_methods_attr {
  $i = 0;
  foreach $method (@method_names) {
    return 1 if ( $method_attr[$i++] eq $_[0] );
  }
  return 0;
}

############################################################################

sub check_members_attr {
  $i = 0;
  foreach $member (@member_names) {
    return 1 if ( $member_attr[$i++] eq $_[0] );
  }
  return 0;
}

############################################################################

sub fads_create_source_file {
  open(OUT, "> $source_file_tmp");
  print(OUT "//\n");
  print(OUT "// File: ${class_root}${src_ext_default}\n");
  print(OUT "//\n");
  print(OUT "//   This file is automatically generated by $script_name.pl.\n");
  print(OUT "//   Do not edit by hand.  Edit $class_name.rootio instead.\n");
  print(OUT "\n");

  print(OUT "#include \"${class_root}${header_ext_default}\"\n\n");

  if ( $assoc_class[0] ne "" ) {
##    print(OUT "// Forward Declarations:\n");
    foreach $aclass (@assoc_class) {
      print(OUT "#include \"${aclass}${header_ext_default}\"\n");
    }
    print(OUT "\n");
  }

  if ( $additional_hdrs_src[0] ne "" ) {
##    print(OUT "// Additional Include:\n");
    foreach $addhdr (@additional_hdrs_src) {
      if ( $addhdr ne "" ) {
        $line = "#include \"$addhdr\"";
        $line =~ s/"</</;
        $line =~ s/>"/>/;
        print(OUT "$line\n");
      }
    }
    print(OUT "\n");
  }

  if ( $root_class == 1 ) {
    print(OUT "ClassImp(${class_root})\n");
    print(OUT "\n");
  }

  if ( $global_impls[0] ne "" ) {
    foreach $impl (@global_impls) {
      print(OUT "$impl\n");
    }
    print(OUT "\n");
  }

  if ( $namespace ne "none" ) {
    $nspace ="${namespace}::";
  } else {
    $nspace = "";
  }

  if ( $constructors[0] ne "" ) {
    $i=0;
    foreach $constructor (@constructors) {
      if ( ! &is_no_implementation($constructor) ) {
        &fads_create_constructor($i++);
      }
    }
  } else {
    &fads_create_constructor(-1);
  }

  if ( $destructors[0] ne "" ) {
    $i=0;
    foreach $destructor (@destructors) {
      if ( ! &is_no_implementation($destructor) ) {
        &fads_create_destructor($i++);
      }
    }
  } else {
    &fads_create_destructor(-1);
  }

  &fads_create_methods;

  print(OUT "// End of ${class_root}${src_ext_default}\n\n");

  close(OUT);

  print "File $source_file created.\n" if ( $verbose && ! $nowrite );
}

############################################################################

sub fads_create_constructor {
  $j = $_[0];
  $k = $j + 1;
  if ( $j >= 0 ) {
    $namarg = $constructors[$j];
    $dname  = "Constructor #" . $k;
    $impl   = $constructor_impl[$j];
  } else {
    $namarg = $class_root . "()";
    $dname  = "Default Constructor";
    $impl   = "";
  }
##  print(OUT "// Implementation of $dname\n");
  print(OUT "${nspace}${class_root}::$namarg\n");

  if ( $impl ne "" ) {
    #
    # explicit implementation
    #
    &fads_create_implementation_body($impl, $num_constructor_impl[$j], "");
  } else {
    #
    # default implementation
    #
    if ( $member_names[0] ne "") {
      print(OUT " : ");
      $j=0;
      foreach $mname (@member_names) {
        $j++;
        print(OUT "${mname}(0)");
        print(OUT ", ") if ( $j != $num_member );
      }
      print(OUT "\n");
    }
    print(OUT "{\n");
    
    foreach $aclass (@assoc_class) {
      print(OUT "  f_${aclass} = new ${nspace}${aclass}();\n") if ($aclass ne "");
    }
    print(OUT "}\n");
  }
  print(OUT "\n");
}

############################################################################

sub count_members {
  $ii = 0;
  foreach $mname (@member_names) {
    $ii++;
  }
  return $ii;
}

############################################################################

sub fads_create_destructor {
  $j = $_[0];
  $k = $j + 1;
  if ( $j >= 0 ) {
    $namarg = $destructors[$j];
    $namarg =~ s/virtual[ ]+//;
    $dname  = "Destructor #" . $k;
    $impl   = $destructor_impl[$j];
  } else {
    $namarg = "~" . $class_root . "()";
    $dname  = "Default Destructor";
    $impl   = "";
  }
##  print(OUT "// Implementation of $dname\n");
  print(OUT "${nspace}${class_root}::$namarg\n");

  if ( $impl ne "" ) {
    #
    # explicit implementation
    #
    &fads_create_implementation_body($impl, $num_destructor_impl[$j], "");
  } else {
    print(OUT "{\n");
    foreach $aclass (@assoc_class) {
      print(OUT "  delete f_${aclass};\n") if ( $aclass ne "" );
    }
    print(OUT "}\n");
  }
  print(OUT "\n");
}

############################################################################

sub fads_create_methods {
  $i=0;
  foreach $method (@method_names) {
    &fads_create_method($i++, $method);
  }
}

############################################################################

sub fads_create_method {
  $namarg = &method_namarg($_[1]);
  ($mname = $namarg) =~ s/\(.*$//;    # chop off arguments
  $rtype = &return_type($_[1]) ; $rtype = "" if ( $rtype eq "none" );
  $rtype =~ s/\<namespace\>/${namespace}::/  if ( $namespace ne "none" );
  $rtype =~ s/\<namespace\>//                if ( $namespace eq "none" );
  $ntype = $namespace . "::"     ; $ntype = "" if ( $namespace eq "none" );

###    print "$_[1]\n";
###    $ss = &is_implicit_inline($_[1]);
###    $tt = &is_inline($_[1]);
###    $uu = &is_virtual($_[1]);
###    print "   implicit_inline = $ss";
###    print "   is_inline = $tt";
###    print "   is_virtual = $uu\n";

  if ( ! &is_no_implementation($_[1]) ) {

##    print(OUT "// Implementation of $mname\n");
    print(OUT "$rtype ${ntype}${class_root}::${namarg}\n");
    $impl = $method_impl[$_[0]];

    if ( $impl ne "" ) {
      #
      # explicit implementation
      #
      &fads_create_implementation_body($impl, $num_method_impl[$_[0]], "");
    } else {
      #
      # default implementation
      #
      print(OUT "{\n");
      if ( $rtype ne "void" ) {
        print(OUT "  // Comment out the next line in actual implementation\n");
        print(OUT "  // return ($rtype) _something_;\n");
      }
      print(OUT "}\n");
    }
    print(OUT "\n");
  }
}

############################################################################

sub fads_print_description {
  @words = quotewords("\n",1,$_[0]);
  foreach $word (@words) {
    print "  // $word\n" if ( $word ne "" );
  }
}

############################################################################

sub fads_create_description {
  $desc = $_[0];
  if ( $desc ne "" ) {
    @words = quotewords("\n",1,$desc);
    foreach $word (@words) {
      print(OUT "      // $word\n") if ( $word ne "" );
    }
  } else {
    print(OUT "      // $_[1]\n");
  }
  print(OUT "\n");
}

############################################################################

sub fads_print_implementation_body {
  @words = quotewords("\n",1,$_[0]);
  $k=0;
  foreach $word (@words) {
    $k++;
    print "  $word\n" if ( $k <= $_[1] );
  }
}

############################################################################

sub fads_create_implementation_body {
  @words = quotewords("\n",1,$_[0]);
  $lspace = $_[2];
  $k=0;
  foreach $word (@words) {
    $k++;
    $word =~ s/^[ ]+$//;   # tidy empty line
    print(OUT "${lspace}$word\n") if ( $k <= $_[1] );
  }
#  print "Debug: fads_create_implementation_body k = $k  n = $_[1]\n";
}

############################################################################

sub fads_create_array_header_file{
  open(OUT, "> $array_header_file_tmp");

  print(OUT "//\n");
  print(OUT "// File: ${class_array}.hh\n");
  print(OUT "//\n");
  print(OUT "//   This file is automatically generated by g4rootio.pl.\n");
  print(OUT "//   Do not edit by hand.  Edit ${class_name}.rootio instead.\n");
  print(OUT "\n");

  $_ = &capitalize(&underscore($class_array));
  print(OUT "#ifndef $_\n");
  print(OUT "#define $_ 1\n");
  print(OUT "\n");

  print(OUT "#include \"TClonesArray.h\"\n");
  print(OUT "#include \"TObject.h\"\n");
  print(OUT "\n");
  print(OUT "class ${class_name};\n");
  print(OUT "\n");
  print(OUT "// Class Description:\n");
  print(OUT "//   Collection class to store each ${class_name}.\n");
  print(OUT "//   Actual store/retrieve operation is done via ${class_io}.\n");
  print(OUT "\n");
  print(OUT "class ${class_array}\n");
  print(OUT " : public TObject\n");
  print(OUT "{\n");
  print(OUT "    public:\n");
  print(OUT "      ${class_array}();\n");
  print(OUT "      ~${class_array}();\n");
  print(OUT "\n");
  print(OUT "    public:\n");
  print(OUT "      bool Add(${class_name}*);\n");
  print(OUT "      int NumEntries() { return fentries; };\n");
  print(OUT "      TClonesArray* GetArray() const {return farray;};\n");
  print(OUT "      void Clear();\n");
  print(OUT "\n");
  print(OUT "    private:\n");
  print(OUT "      Int_t fentries;          // number of entries;\n");
  print(OUT "      TClonesArray* farray;    //-> \"${class_root}\";\n");
  print(OUT "      static TClonesArray* fgarray;\n");
  print(OUT "\n");
  print(OUT "    ClassDef (${class_array},1)\n");
  print(OUT "\n");
  print(OUT "}; // End of class ${class_array}\n");
  print(OUT "\n");
  print(OUT "#endif\n");

  close(OUT);
  print "File $array_header_file created.\n" if ( $verbose && ! $nowrite );
}

############################################################################

sub fads_create_array_source_file{
  open(OUT, "> $array_source_file_tmp");

  print(OUT "//\n");
  print(OUT "// File: ${class_array}.cc\n");
  print(OUT "//\n");
  print(OUT "//   This file is automatically generated by g4rootio.pl.\n");
  print(OUT "//   Do not edit by hand.  Edit ${class_name}.rootio instead.\n");
  print(OUT "\n");

  print(OUT "#include \"${class_array}.hh\"\n");
  print(OUT "#include \"${class_name}.hh\"\n");
  print(OUT "#include \"${class_root}.hh\"\n");
  print(OUT "\n");
  print(OUT "ClassImp(${class_array})\n");
  print(OUT "\n");
  print(OUT "TClonesArray* ${class_array}::fgarray = 0;\n");
  print(OUT "\n");
  print(OUT "${class_array}::${class_array}()\n");
  print(OUT "{\n");
  print(OUT "  if (!fgarray) fgarray = new TClonesArray(\"${class_root}\",1000);\n");
  print(OUT "  farray = fgarray;\n");
  print(OUT "  fentries = 0;\n");
  print(OUT "}\n");
  print(OUT "\n");
  print(OUT "${class_array}::~${class_array}()\n");
  print(OUT "{\n");
  print(OUT "  Clear();\n");
  print(OUT "}\n");
  print(OUT "\n");
  print(OUT "bool ${class_array}::Add(${class_name}* a)\n");
  print(OUT "{\n");
  print(OUT "  TClonesArray &ary = *farray;\n");
  print(OUT "  ${class_root}* h = new(ary[fentries++]) ${class_root}(a);\n");
  print(OUT "  if (h==0) return false;\n");
  print(OUT "  return true;\n");
  print(OUT "}\n");
  print(OUT "\n");
  print(OUT "void ${class_array}::Clear()\n");
  print(OUT "{\n");
  print(OUT "  fentries = 0;\n");
  print(OUT "  farray->Delete();\n");
  print(OUT "}\n");
  print(OUT "\n");
  print(OUT "// End of ${class_array}.cc\n");

  close(OUT);
  print "File $array_source_file created.\n" if ( $verbose && ! $nowrite );
}

############################################################################

sub fads_create_io_header_file{
  open(OUT, "> $io_header_file_tmp");

  print(OUT "//\n");
  print(OUT "// File: ${class_io}.hh\n");
  print(OUT "//\n");
  print(OUT "//   This file is automatically generated by g4rootio.pl.\n");
  print(OUT "//   Do not edit by hand.  Edit ${class_name}.rootio instead.\n");
  print(OUT "\n");

  $_ = &capitalize(&underscore($class_io));
  print(OUT "#ifndef $_\n");
  print(OUT "#define $_ 1\n");
  print(OUT "\n");

  print(OUT "#include <string>\n");
  print(OUT "#include \"G4HitRootIO.hh\"\n");
  print(OUT "#include \"${class_name}.hh\"\n");
  print(OUT "#include \"${class_array}.hh\"\n");
  print(OUT "#include \"G4RootTransManager.hh\"\n");
  print(OUT "#include \"${array_io_base}.hh\"\n");
  print(OUT "\n");
  print(OUT "class ${class_array};\n");
  print(OUT "\n");
  print(OUT "class ${class_io}\n");
  print(OUT " : public ${array_io_base}\n");
  print(OUT "{\n");
  print(OUT "    public:\n");
  print(OUT "      ${class_io}(G4std::string detName, G4std::string colName);\n");
  print(OUT "      ~${class_io}();\n");
  print(OUT "\n");
  print(OUT "    public:\n");
  print(OUT "      bool Store(const ${collection_base_class}* hc);\n");
  print(OUT "      bool Retrieve( ${collection_base_class}*& hc);\n");
  print(OUT "\n");
  print(OUT "    private:\n");
  print(OUT "      G4std::string       f_branchName;\n");
  print(OUT "      G4std::string       f_branchDesc;\n");
  print(OUT "      G4RootTransManager* f_transMan;\n");
  print(OUT "      ${class_array}* f_hc;\n");
  print(OUT "      int f_nev;\n");
  print(OUT "\n");
  print(OUT "}; // End of class ${class_io}\n");
  print(OUT "\n");
  print(OUT "#endif\n");

  close(OUT);
  print "File $io_header_file created.\n" if ( $verbose && ! $nowrite );
}

############################################################################

sub fads_create_io_source_file{
  open(OUT, "> $io_source_file_tmp");

  print(OUT "//\n");
  print(OUT "// File: ${class_io}.cc\n");
  print(OUT "//\n");
  print(OUT "//   This file is automatically generated by g4rootio.pl.\n");
  print(OUT "//   Do not edit by hand.  Edit ${class_name}.rootio instead.\n");
  print(OUT "\n");
  print(OUT "#include \"${class_io}.hh\"\n");
  print(OUT "#include \"${collection_class}.hh\"\n") if ( $collection_class_header );
  print(OUT "#include \"${catalog}.hh\"\n");
  print(OUT "#include \"G4PersistencyManager.hh\"\n");
  print(OUT "#include \"G4RootIOManager.hh\"\n");
  print(OUT "#include \"${class_root}.hh\"\n");
  print(OUT "#include \"TFile.h\"\n");
  print(OUT "#include \"TTree.h\"\n");
  print(OUT "#include \"TBranch.h\"\n");
  print(OUT "\n");
  print(OUT "static ${catalog}<${class_io}> ${entry_object}(\"${sdet_name}\");\n");
  print(OUT "\n");
  print(OUT "${class_io}::${class_io}(G4std::string detName, G4std::string colName)\n");
  print(OUT " : ${array_io_base}(detName, colName), f_nev(0)\n");
  print(OUT "{\n");
  print(OUT "  f_branchName = colName + \".\";\n");
  print(OUT "\n");
  print(OUT "  G4RootIOManager* pm = (G4RootIOManager*)\n");
  print(OUT "     G4PersistencyCenter::GetPersistencyCenter()->CurrentPersistencyManager();\n");
  print(OUT "  assert(pm!=0);\n");
  print(OUT "  f_transMan = (G4RootTransManager*) pm->TransactionManager();\n");
  print(OUT "  assert(f_transMan!=0);\n");
  print(OUT "\n");
  print(OUT "  // create a placeholder for new persistent array\n");
  print(OUT "  f_hc = new ${class_array}();\n");
  print(OUT "  assert(f_hc!=0);\n");
  print(OUT "\n");
  print(OUT "  // Tell ROOT the class name to be stored in this branch\n");
  print(OUT "  f_branchDesc = \"${class_array}\";\n");
  print(OUT "}\n");
  print(OUT "\n");
  print(OUT "${class_io}::~${class_io}()\n");
  print(OUT "{\n");
  print(OUT "  delete f_hc;\n");
  print(OUT "}\n");
  print(OUT "\n");
  print(OUT "bool ${class_io}::Store(const ${collection_base_class}* hc)\n");
  print(OUT "{\n");
  print(OUT "  // cast the transient collection for this detector type.\n");
  print(OUT "  ${collection_class}* ahc = (${collection_class}*) hc;\n");
  print(OUT "  if ( ahc->$collection_class_method_entries <= 0 ) {\n");
  print(OUT "    return true;\n");
  print(OUT "  }\n");
  print(OUT "\n");
  print(OUT "  f_transMan->LastWriteFile()->cd();\n");
  print(OUT "\n");
  print(OUT "  int bufsize = f_transMan->BufferSize();\n");
  print(OUT "  int split = f_transMan->SplitLevel();\n");
  print(OUT "\n");
  print(OUT "  // Set the new branch style\n");
  print(OUT "  TTree::SetBranchStyle(1);\n");
  print(OUT "\n");
  print(OUT "  TTree* tree = f_transMan->LastWriteTree();\n");
  print(OUT "\n");
  print(OUT "  // Get the branch if already exists\n");
  print(OUT "  TBranch* branch = tree->GetBranch(f_branchName.c_str());\n");
  print(OUT "\n");
  print(OUT "  // Set the branch with \"collection_name.class_name\"\n");
  print(OUT "  if ( branch == 0 ) {\n");
  print(OUT "    branch = tree->Branch(f_branchName.c_str(), f_branchDesc.c_str(),\n");
  print(OUT "                          &f_hc, bufsize, split);\n");
  print(OUT "    if ( m_verbose > 2 ) {\n");
  print(OUT "      G4cout << \"${class_io}: Branch \\\"\" << f_branchName\n");
  print(OUT "                << \"\\\" initialized with bufsize = \" << bufsize\n");
  print(OUT "                << \", split level = \" << split << \".\" << G4endl;\n");
  print(OUT "    }\n");
  print(OUT "  } else {\n");
  print(OUT "    branch->SetAddress(&f_hc);\n");
  print(OUT "  }\n");
  print(OUT "  assert(branch != 0);\n");
  print(OUT "  branch->SetAutoDelete(kFALSE);\n");
  print(OUT "\n");
  print(OUT "  if ( m_verbose > 2 ) {\n");
  print(OUT "    G4cout << \"${class_io}: Storing ${class_name} collectoin \"\n");
  print(OUT "              << \"for \\\"\" << ahc->GetSDname() << \"\\\", \"\n");
  print(OUT "              << \"collection \\\"\" << ahc->GetName() << \"\\\".\" << G4endl\n");
  print(OUT "              << \"  # of entries: \" << ahc->$collection_class_method_entries << \".\" << G4endl;\n");
  print(OUT "  }\n");
  print(OUT "\n");
  print(OUT "  // reset the persistent array\n");
  print(OUT "  f_hc->Clear();\n");
  print(OUT "\n");
  print(OUT "  // loop for the transient collections and copy the values to\n");
  print(OUT "  // persistent collection\n");
  print(OUT "  for ( int i = 0; i < ahc->$collection_class_method_entries; i++ ) {\n");
  print(OUT "    ${class_name}* a = (*ahc)$collection_class_method_get;\n");
  print(OUT "    if ( a != 0 ) {\n");
  print(OUT "      if ( m_verbose > 1 )\n");
  print(OUT "        G4cout << \"${class_io}: Storing ${class_name} #\"\n");
  print(OUT "                  << i << \".\" << G4endl;\n");
  print(OUT "       f_hc->Add(a);\n");
  print(OUT "    }\n");
  print(OUT "  }\n");
  print(OUT "\n");
  print(OUT "  tree->Fill();\n");
  print(OUT "\n");
  print(OUT "  return true;\n");
  print(OUT "}\n");
  print(OUT "\n");
  print(OUT "bool ${class_io}::Retrieve( ${collection_base_class}*& hc)\n");
  print(OUT "{\n");
  print(OUT "  hc = 0;\n");
  print(OUT "\n");
  print(OUT "  // Get the ROOT tree \"Geant4\"\n");
  print(OUT "  TTree* tree = f_transMan->LastReadTree();\n");
  print(OUT "  if ( tree == 0 ) {\n");
  print(OUT "    G4cerr << \"${class_io}::Retrieve: error getting the ROOT file tree \\\"Geant4\\\".\" << G4endl;\n");
  print(OUT "    return false;\n");
  print(OUT "  }\n");
  print(OUT "\n");
  print(OUT "  // Set the branch with \"collection_name.class_name\"\n");
  print(OUT "  TBranch* branch = tree->GetBranch(f_branchName.c_str());\n");
  print(OUT "  if ( branch == 0 ) {\n");
  print(OUT "    G4cerr << \"${class_io}::Retrieve: error getting the ROOT file branch \\\"\" << f_branchName << \"\\\".\" << G4endl;\n");
  print(OUT "    return false;\n");
  print(OUT "  }\n");
  print(OUT "  branch->SetAddress(&f_hc);\n");
  print(OUT "\n");
  print(OUT "  if ( m_verbose > 2 ) {\n");
  print(OUT "    G4cout << \"${class_io}: Retrieving ${class_name} collectoin \"\n");
  print(OUT "              << \"for \\\"\" << f_detName << \"\\\", \"\n");
  print(OUT "              << \"collection \\\"\" << f_colName << \"\\\".\" << G4endl;\n");
  print(OUT "  }\n");
  print(OUT "\n");
  print(OUT "  // Read an n-th \"event\" from the tree\n");
  print(OUT "  Int_t nb = tree->GetEntry(f_nev++);\n");
  print(OUT "  if ( nb == 0 ) return false;\n");
  print(OUT "\n");
  print(OUT "  if ( m_verbose > 4 )\n");
  print(OUT "    G4cout << \"${class_io}: Retrieved \" << nb\n");
  print(OUT "              << \" bytes for event \" << f_nev-1 << \".\" << G4endl;\n");
  print(OUT "\n");
  print(OUT "  // create a new transient collection\n");
  print(OUT "  ${collection_class}* ahc =\n");
  print(OUT "               new ${collection_class_signature};\n");
  print(OUT "\n");
  print(OUT "  // Initialize the iterator\n");
  print(OUT "  TClonesArray* ary = f_hc->GetArray();\n");
  print(OUT "  TIter NextEntry(ary);\n");
  print(OUT "\n");
  print(OUT "  // Loop over each entry\n");
  print(OUT "  int nh = 0;\n");
  print(OUT "  while ( ${class_root}* h = (${class_root}*)NextEntry() ) {\n");
  print(OUT "    nh++;\n");
  print(OUT "    ${class_name}* a = h->${make_transient}();\n");
  print(OUT "    ahc->$collection_class_method_add;\n");
  print(OUT "  }\n");
  print(OUT "\n");
  print(OUT "  // reset the persistent array\n");
  print(OUT "  f_hc->Clear();\n");
  print(OUT "\n");
  print(OUT "  // Check the number of stored entries\n");
  print(OUT "  return true;\n");
  print(OUT "}\n");
  print(OUT "\n");
  print(OUT "// End of ${class_io}.cc\n");

  close(OUT);
  print "File $io_source_file created.\n" if ( $verbose && ! $nowrite );
}

############################################################################

sub is_virtual {
  @words = quotewords(" ", 1, $_[0]);
  $result = 0;
  foreach $word (@words) {
    $result = 1 if ( $word eq "virtual");
  }
  return $result;
}

############################################################################

sub is_static {
  @words = quotewords(" ", 1, $_[0]);
  $result = 0;
  foreach $word (@words) {
    $result = 1 if ( $word eq "static");
  }
  return $result;
}

############################################################################

sub is_inline {
  @words = quotewords(" ", 1, $_[0]);
  $result = 0;
  foreach $word (@words) {
    $result = 1 if ( $word eq "inline");
  }
  return $result;
}

############################################################################

sub is_implicit_inline {
  $result = 0;
  $result = 1 if ( $_[0] =~ /{.*}/ );
  return $result;
}

############################################################################

sub is_pure {
  $_ = $_[0];
  s/(^.*)(=[ ]*0[ ]*$)/$2/;
  $result = 0;
  $result = 1 if ( /=[ ]*0[ ]*$/ );
  return $result;
}

############################################################################

sub is_no_implementation {
#  if ($verbose ) {
#    $st1 = &is_implicit_inline($_[0]);
#    $st2 = &is_inline($_[0]);
#    $st3 = &is_virtual($_[0]);
#    print "method $_[0]  st = $st1  $st2  $st3\n";
#  }
  $result = 0;
  $result = 1 if ( &is_implicit_inline($_[0]) || &is_inline($_[0]) ||
                   &is_pure($_[0]) );
#                   &is_virtual($_[0]) );
#    print "method $_[0]  result = $result\n" if ($verbose);
  return $result;
}

############################################################################

sub is_endline {
  $_ = $_[0] if ( $_[0] ne "" );
  return(1) if ( /^\.\.$/ );
  return(0);
}

############################################################################

sub is_description {
  $ll = $_[0];
  $ll =~ s/^#.*$//;   # omit comments
  $ll =~ s/^[ ]+//;   # chop off leading spaces
  if ( $ll =~ /^\/\// ) {
    print "is_description: return 1\n";
    return(1) if ( $ll =~ /^\/\// );
  }
  return(0);
}

############################################################################

sub is_no_source {
  $result = 1;
  if ( $template ne "" ) {
    return(1);
  }
  if ( $constructors[0] ne "" ) {
    foreach $constructor (@constructors) {
      return(0) if ( ! &is_no_implementation($constructor) );
    }
  } else {
    return(0);
  }
  if ( $destructors[0] ne "" ) {
    foreach $destructor (@destructors) {
      return(0) if ( ! &is_no_implementation($destructor) );
    }
  } else {
    return(0);
  }
  if ( $method_names[0] ne "" ) {
    foreach $method (@method_names) {
      return(0) if ( ! &is_no_implementation($method) );
    }
  }
  return $result;
}

############################################################################

sub return_type {
  @words = quotewords(" ", 1, $_[0]);
  $n=0;
  foreach $word (@words) {
    $n++ if ( $word eq "virtual");
    $n++ if ( $word eq "static");
    $n++ if ( $word eq "inline");
  }
  return $words[$n];
}

############################################################################

sub method_namarg {
  @words = quotewords(" ", 1, $_[0]);
  $n=0;
  foreach $word (@words) {
    shift @words if ( $word eq "virtual");
    shift @words if ( $word eq "static");
    shift @words if ( $word eq "inline");
    $n++;
  }
  shift @words;
  $result = "";
  $ii=0;
  foreach $word (@words) {
    $result .= $word;
    $result .= " " if ( $ii++ < $n );
  }
  $result =~ s/[ ]+$//;  # chop off the trailing branks
  return $result;
}

############################################################################

sub check_dir {
  while ( ! -d $_[0] ) {
    if ( ! $force ) {
      print "Directory $_[0] does not exist.  Create it? (yes/no) : ";
      <IN>;
      chop;
      print;
      if ( $_ eq "y" || $_ eq "yes" ) {
        system("mkdir -p $_[0]");
      } elsif ( $_ eq "n" || $_ eq "no" ) {
        print "exiting...\n";
        exit 1;
      } elsif ( $_ ne "" ) {
        print "\nError - answer yes or no!\n";
      }
    } else {
      print "\n${script_name}: Error - directory $_[0] does not exist.  Stop.\n";
      exit 1;
    }
  }
}
  
############################################################################

sub check_file {
  if ( ! $force ) {
    while ( -f $_[0] ) {
      print "File $_[0] already exists - discard old file? (yes/no) : ";
      <IN>;
      chop;
      if ( $_ eq "n" || $_ eq "no" ) {
        print "\nSpecify alternative header file name: ";
        $_[0] = <IN>;
        chop($_[0]);
        print "\nNew header file name: $_[0]\n";
      } elsif ( $_ eq "y" || $_ eq "yes" ) {
        system("rm -f $_[0]");
      } elsif ( $_ ne "" ) {
        print "\nError - answer yes or no!\n";
      }
    }
  }
}

############################################################################

sub diff_file {
  if ( -f $_[1] ) {
    if ( -f $_[0] ) {
      system ("egrep -v '^// \\\$Id\: | tag \\\$Name\: ' $_[0] > $_[0].tmp");
      system ("egrep -v '^// \\\$Id\: | tag \\\$Name\: ' $_[1] > $_[1].tmp");
      system("diff $_[1].tmp $_[0].tmp >& /dev/null");
      if ( $? ) {
##      if ( compare($_[0].tmp, $_[1].tmp) == 0 ) {
        print  "diff $_[1].tmp $_[0].tmp\n";
        system("diff $_[1].tmp $_[0].tmp");
      }
      system ("rm -rf $_[0].tmp $_[1].tmp");
    }
  } else {
    print "Original file $_[1] does not exist.\n";
  }
}

############################################################################

sub install_file {
  if ( -f $_[0] ) {
    if ( -f $_[1] ) {
      system ("egrep -v '^// \\\$Id\: | tag \\\$Name\: ' $_[0] > $_[0].tmp");
      system ("egrep -v '^// \\\$Id\: | tag \\\$Name\: ' $_[1] > $_[1].tmp");
      system("diff $_[1].tmp $_[0].tmp >& /dev/null");
      if ( $? ) {
##      if ( compare($_[0].tmp, $_[1].tmp) != 0 ) {
        print "Replacing $_[1] ... \n";
        system ("mv -f $_[0] $_[1]");
        exit 1 if $?;
##      } else {
##        print "File $_[1] and $_[0] are same.\n";
      }
      system ("rm -rf $_[0].tmp $_[1].tmp");
    } else {
        print "Installing $_[1] ... \n";
        system ("mv -f $_[0] $_[1]");
        exit 1 if $?;
    }
  } else {
    print "File $_[0] does not exist!\n";
    exit 1;
  }
}

############################################################################
