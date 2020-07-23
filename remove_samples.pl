use strict ; 
use warnings ; 

my %samples ; 
open IN, "<samples_in_latest_tree.txt" ;
while (<IN>) { 
	chomp ; 
	$samples{$_} ++ ; 
}
close IN ; 
my @keep ;
open IN, '<', $ARGV[0] ;
while (<IN>) { 
	chomp ; 
	my @split = split ( /\t/, $_ ) ; 
	if ( $split[0] =~ m/\#CHROM/ ) {
		foreach ( 9..$#split ) {
                  
			if ( exists( $samples{$split[$_]} ) ) {
				push @keep, $_ ; 
			}
		}
		last ; 
	}	
}
close IN ; 

open IN, '<', $ARGV[0] ;
while (<IN>) {
        chomp ;
	if ( $_ =~ m/\#\#/ ) { 
		print $_, "\n" ; 
	}
	else {
	        my @split = split ( /\t/, $_ ) ;
		print $split[0] ;
                foreach ( 1..8,@keep ) { 
			print "\t", $split[$_] ;
                }
		print "\n" ;
        }
}
close IN ; 
