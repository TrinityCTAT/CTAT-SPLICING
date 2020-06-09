#!/usr/bin/env perl

use strict;
use warnings;

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use TiedHash;


my $usage = <<__EOUSAGE__;

#####################################################################
#
# --cancer_introns_tsv <string>   cancer introns tsv file
#
# --ctat_genome_lib <string>      /path/to/ctat_genome_lib_build_dir
#
#####################################################################

__EOUSAGE__

    ;


my $help_flag;

my $cancer_introns_tsv;
my $ctat_genome_lib;



&GetOptions ( 'h' => \$help_flag,

              'cancer_introns_tsv=s' => \$cancer_introns_tsv,
              'ctat_genome_lib=s' => \$ctat_genome_lib,
    );

if ($help_flag) { die $usage; }


unless ($cancer_introns_tsv && $ctat_genome_lib) {
    die $usage;
}


main: {

    
    my $cancer_introns_db_dir = "$ctat_genome_lib/cancer_splicing_lib";
    if (! -d $cancer_introns_db_dir) {
        mkdir($cancer_introns_db_dir) or die "Error, could not mkdir $cancer_introns_db_dir";
    }
    
    my $db_idx_file = "$cancer_introns_db_dir/cancer_splicing.idx";
        
    my $idx = new TiedHash( { create => $db_idx_file } );
    
    ## store placeholder to ensure it's working in applications that leverage it.
    $idx->store_key_value("chr:ABC-DEF", "__placeholder_testval__");
    
    my $fh;
    if ($cancer_introns_tsv =~ /\.gz$/) {
        open($fh, "gunzip -c $cancer_introns_tsv | ") or die "Error, cannot open file $cancer_introns_tsv via gunzip ";
    }
    else {
        open($fh, $cancer_introns_tsv) or die "Error, cannot open file: $cancer_introns_tsv";
    }
    
    my $column_header_line = <$fh>;
    chomp $column_header_line;
    $idx->store_key_value("column_headers", $column_header_line);
    
    while (<$fh>) {
        chomp;
        my ($intron, $annot_text) = split(/\t/, $_, 2);
        
        $idx->store_key_value($intron, $annot_text);
        
    }
    close $fh;
    
    print STDERR "\n\nDone building cancer splicing db: $db_idx_file\n";
    
    exit(0);
}
           
