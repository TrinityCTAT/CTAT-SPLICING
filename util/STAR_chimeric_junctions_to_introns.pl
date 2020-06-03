#!/usr/bin/env perl

# contributed by Brian Haas, Broad Institute, 2015

use strict;
use warnings;
use Carp;
use Cwd;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use File::Basename;
use Data::Dumper;
use ChimericCigarParser;


use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

## Options
my $chimeric_junction_file;
my $help_flag;


my $usage = <<__EOUSAGE__;

###################################################################################
#
#  Required:
#
#    --chimeric_junction|J <string>     Chimeric.out.junction file
#
###################################################################################


__EOUSAGE__

    ;


my $MIN_INTRON_LEN = 20;

my $DEBUG;

&GetOptions ( 'h' => \$help_flag,
              
              'chimeric_junction|J=s' => \$chimeric_junction_file,
              
              'd' => \$DEBUG,
    );


if ($help_flag) {
    die $usage;
}
unless ($chimeric_junction_file) {
    die $usage;
}


main: {

    &map_chimeric_reads_to_introns($chimeric_junction_file);
    
    exit(0);

    
}


####
sub map_chimeric_reads_to_introns {
    my ($junctions_file) = @_;
    
    print STDERR "-mapping reads to genes\n";

    
    if ($junctions_file =~ /\.gz$/) {
        $junctions_file = "gunzip -c $junctions_file | ";
    }


    my %intron_counter;

    my $start_time = time();
    my $counter = 0;
    open (my $fh, $junctions_file) or die "Error, cannot open file $junctions_file";
    while (<$fh>) {
        if (/^\#/) { next; }
        if (/^chr_donorA\tbrkpt/) { 
            # header line
            next;
        }
        $counter++;
        if ($counter % 100 == 0) {
            my $time = time();
            my $seconds = $time - $start_time;
            if ($seconds && $counter % 10000 == 0) {
                my $rate = sprintf("%.2f", $counter / $seconds * 60);
                print STDERR "\r[$counter], rate=$rate/min ";
            }
        }
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        
        # from star doc:
        #The rst 9 columns give information about the chimeric junction:

        #The format of this le is as follows. Every line contains one chimerically aligned read, e.g.:
        #chr22 23632601 + chr9 133729450 + 1 0 0 SINATRA-0006:3:3:6387:56650 23632554 47M29S 133729451 47S29M40p76M
        #The first 9 columns give information about the chimeric junction:

        #column 1: chromosome of the donor
        #column 2: rst base of the intron of the donor (1-based)
        #column 3: strand of the donor
        #column 4: chromosome of the acceptor
        #column 5: rst base of the intron of the acceptor (1-based)
        #column 6: strand of the acceptor
        #column 7: junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC
        #column 8: repeat length to the left of the junction
        #column 9: repeat length to the right of the junction
        #Columns 10-14 describe the alignments of the two chimeric segments, it is SAM like. Alignments are given with respect to the (+) strand
        #column 10: read name
        #column 11: rst base of the rst segment (on the + strand)
        #column 12: CIGAR of the rst segment
        #column 13: rst base of the second segment
        #column 14: CIGAR of the second segment

        my $junction_type = $x[6];
                
        my $read_name = $x[9];
        
        my ($chrA, $coordA, $orientA) = ($x[0], $x[1], $x[2]);
        $coordA = ($orientA eq '+') ? --$coordA : ++$coordA;
        
        my ($rst_A, $cigar_A) = ($x[10], $x[11]);
        

        my ($chrB, $coordB, $orientB) = ($x[3], $x[4], $x[5]);
        $coordB = ($orientB eq '+') ? ++$coordB : --$coordB;
        
        
        unless ($chrA eq $chrB) { next; } # want introns for long introns and readthru splices only.
        
        unless ($orientA eq $orientB) { next; } # want consistent orientations.

        my ($rst_B, $cigar_B) = ($x[12], $x[13]);

        my ($genome_coords_A_aref, $read_coords_A_aref) = &get_genome_coords_via_cigar($rst_A, $cigar_A);
        my ($genome_coords_B_aref, $read_coords_B_aref) = &get_genome_coords_via_cigar($rst_B, $cigar_B);
        
        #print "$chrA\t$genome_coords_A_aref->[0]->[0]\n";
        
        
        
        my @genomic_coords_A = @$genome_coords_A_aref;
        my @genomic_coords_B = @$genome_coords_B_aref;

        my @genomic_coords;
        if ($orientA eq '+') {
            @genomic_coords = (@genomic_coords_A, @genomic_coords_B);
        }
        else {
            @genomic_coords = (@genomic_coords_B, @genomic_coords_A);
        }

               
        for (my $i = 0; $i < $#genomic_coords; $i++) {
            my $intron_lend = $genomic_coords[$i]->[1] += 1;
            my $intron_rend = $genomic_coords[$i+1]->[0] - 1;

            if ($intron_rend - $intron_lend > $MIN_INTRON_LEN) {
                
                $intron_counter{"$chrA:$intron_lend-$intron_rend"}++;
            }
        }
        
    }
    
    foreach my $intron (keys %intron_counter) {
        my $count = $intron_counter{$intron};
        
        print join("\t", $intron, $count) . "\n";
    }
    
    
    return;
    
}

