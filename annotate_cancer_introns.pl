#!/usr/bin/env perl

use strict;
use warnings;

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use TiedHash;


my $ctat_genome_lib = $ENV{CTAT_GENOME_LIB};


my $usage = <<__EOUSAGE__;

#####################################################################
#
# --introns_file <string>          cancer introns tsv file (requires a columns header line)
#
# --ctat_genome_lib <string>      /path/to/ctat_genome_lib_build_dir (default: $ctat_genome_lib)
#
# optional:
#
# --intron_col <int>              tab-delim column index for intron (default: 0)
#
#####################################################################

__EOUSAGE__

    ;


my $help_flag;

my $introns_file;
my $intron_col = 0;


&GetOptions ( 'h' => \$help_flag,

              'introns_file=s' => \$introns_file,
              'ctat_genome_lib=s' => \$ctat_genome_lib,
              
              'intron_col=i' => \$intron_col,
    );


if ($help_flag) { die $usage; }


unless ($introns_file && $ctat_genome_lib) {
    die $usage;
}


main: {

    
    my $cancer_introns_db_dir = "$ctat_genome_lib/cancer_splicing_lib";
    my $db_idx_file = "$cancer_introns_db_dir/cancer_splicing.idx";
    
    unless (-s $db_idx_file) {
        die "Error, cannot locate resource file: $db_idx_file";
    }
    
    my $idx = new TiedHash( { 'use' => $db_idx_file } );
    

    ## store placeholder to ensure it's working in applications that leverage it.
    unless($idx->get_value("chr:ABC-DEF") eq "__placeholder_testval__") {
        die "Error, $db_idx_file doesn't appear useable and must be rebuilt";
    }



    
    my $found = 0;
    open(my $fh, $introns_file) or die "Error, cannot open file: $introns_file";

    # print header:
    my $header = <$fh>;
    chomp $header;
    my $header_add = $idx->get_value("column_headers");
    print join("\t", $header, $header_add) . "\n";
    
    while (<$fh>) {
        chomp;
        my $input_line = $_;
        my @x = split(/\t/);
        my $intron = $x[$intron_col];
                
        my $intron_annot = $idx->get_value($intron);
        if ($intron_annot) {
            print join("\t", $input_line, $intron_annot) . "\n";
            $found += 1;
        }
    }
    close $fh;
    
    if ($found) {
        print STDERR "-$introns_file: identified $found cancer introns\n";
    }
    else {
        print STDERR "-$introns_file: no cancer introns identified.\n";
    }
    
    exit(0);
}
           
