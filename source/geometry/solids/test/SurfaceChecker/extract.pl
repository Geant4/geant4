#!/usr/bin/perl -w 

use strict ;

use File::Basename ;

# variables 
my ($i,$j)  ;              # loop
my $events ;         # number of entries 

my @ret ;  ;   # variable handling (for i/o ntuple)

my $mcount = 0;
my $txtfile ;
my $hbkfile ;
my $file ;

my @a ;
my @surface = () ;  # surface point  XX (generated)
my @vertex = (0,0,0) ;   # vertex position P (generated)
my @momentum = () ; # momentum p (generated)
my @intersection = () ;  # reconstructed intersection with solid
my $delta = 0 ;       # distance between surface point XX and intersection
my @bestintersection = () ;  # the best intersection (there are more than 1) 
my $bestdelta = 0 ;
my $phi = 0 ;  # two parameters describing the surface point 
my $u  = 0;

my $distance ;  # distance between surface point and true vertex
my $theta ;     # angle between momentum p and surface normal at XX (true)

# --------------------------------------------------------


# now start a loop over different approximations ..

my @filelist = glob "data/*.data" ;

foreach $file ( @filelist ) {

    next unless ( -s $file ) ; 
    
    $txtfile = "hbk/" . basename($file,".data") . ".txt" ; 
    $hbkfile = "hbk/" . basename($file,".data") . ".hbk" ;
 
    if ( -e $hbkfile )
    {
	print "merge of $hbkfile already done...\n" ;
	next ;
    }

    print "merge $file to temporary $txtfile\n" ;

    open(TXT,">$txtfile") || die "cannot open file $txtfile:$!\n";

    open(FH,$file) || die "cannot open file $file:$!\n";
    $mcount = 0 ;
    while ( <FH> ) {

	next unless /^(Event|Surface|Vertex|Momentum|EndOfEvent|Intersection|Distance|Angle)/ ; 

	s/[(),]/ /g ;            # substitute brackets {} and ,
	@a = split ;             # split the line

	if ( $a[0] eq "Event" ) { 

	    $mcount ++ ;    # event counter

	    $bestdelta = 1e99 ;   # reset the data for every event
	    @surface = ( 1e99,1e99,1e99) ; 
	    @vertex = (1e99,1e99,1e99) ;
	    @momentum = (1e99,1e99,1e99) ;
	    @intersection = (1e99,1e99,1e99) ;
	    $distance = 1e99 ;
	    $theta = 1e99 ;

	    if ( $mcount % 10000 == 0 ) {
		print "event nr. $mcount\n" ;
	    }
	}
	elsif ( $a[0] eq "Surface" ) {
	    @surface = (  $a[1], $a[2], $a[3] )  ;
	}
	elsif ( $a[0] eq "Vertex" ) {
	    @vertex = ( $a[1], $a[2], $a[3] )  ;
	    $u = $a[4] ; $phi = $a[5] ;
	}
	elsif ( $a[0] eq "Momentum" ) {
	    @momentum = ( $a[1], $a[2], $a[3] )  ;
	}
	elsif ( $a[0] eq "Distance" ) {
	    $distance = $a[1] ;
	}
	elsif ( $a[0] eq "Angle" ) {
	    $theta = $a[1] ;
	}
	elsif ( $a[0] eq "Intersection" ) {
	    @intersection = ( $a[1], $a[2], $a[3] )  ;
	    $delta = $a[4] ;
	    if ( $delta < $bestdelta ) { 
		$bestdelta = $delta ;
		@bestintersection = @intersection ;
	    }
	}
	elsif ( $a[0] eq "EndOfEvent" ) {
	    @ret = ( $phi, $u, @surface , @vertex, 
		     $distance, @momentum, @bestintersection, 
		     $bestdelta, $theta)  ;

	    print TXT "@ret\n" ;
	}



    }

    print "$mcount events merged to file $txtfile\n" ;

    close FH ;
    close TXT ;

    print "start converting ascii $txtfile to hbook file $hbkfile...\n" ;

    system("export arg1=$txtfile; export arg2=$hbkfile; export arg3=$mcount ; paw -w 0 -b ntuple") && die "cannot create the ntuple in $hbkfile:$!\n" ;

    # delete the temporary ascii file. Comment out this line if needed.
    system("rm $txtfile") && die "cannot remove $txtfile.$!\n" ;

    print "\n\n" ;

}

