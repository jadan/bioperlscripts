#!/bin/perl -w #!/usr/bin/perl
#Base script modified from https://www.biostars.org/p/46281/
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

$| = 1;

print "Usage: perl extractintergenic.pl fastafile.fasta gff3file.gff3 nameofoutput";

my $outfile_intergenic = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].intergenic.fasta" );

# First, index the genome
my $file_fasta = $ARGV[0];
my $db = Bio::DB::Fasta->new($file_fasta);
print ("Genome fasta parsed\n");

# Second, parse the GFF3
my $last_gene_end = 0;
my $last_gene_frame = 'o';
my $last_gene_id = 'o';

# counter: [0]= contains 1, [1] = contains 2, [2] = doesn't contain  
my @counter = (0,0,0);
my $total_sequences_extracted = 0;
my $total_length = 0;

open GFF, "<$ARGV[1]" or die $!;

while ( my $line = <GFF> ) {
    chomp $line;
    my @array = split( "\t", $line );
    my $type = $array[2];

    #Skip if GFF line is a comment, the genes are in a scaffold, or if the last gene was the last one in a chromosome.
    if (($array[0] =~ m/#/) or ($array[0] =~ m/scaffold/)) {
        next;
    }
    if($last_gene_id ne $array[0]){
        next;
    }

    if (($type eq 'gene') and ( $last_gene_frame ne 'o' ) ) {
        $total_sequences_extracted += 1;

        my $current_frame = $array[6];
        my $gene_name = "prueba".$total_sequences_extracted;
        my $gene_start = $array[3];
        my $gene_end = $array[4];

        # The intergenic sequence.
        my $intergenic_start = $last_gene_end + 1;
        my $intergenic_end = $gene_start - 1;

        my $intergenic_seq = $db->seq( $array[0], $intergenic_start, $intergenic_end);
        
        my $output_intergenic = Bio::Seq->new(
            -seq        => $intergenic_seq,
            -id         => $gene_name."_intergenic",
            -display_id => $gene_name."_intergenic",
            -alphabet   => 'dna',
        );

        if($intergenic_end - $intergenic_start < 10000){
            print($gene_name ." " . $output_intergenic->length() . "\n")
        }

        $total_length += $output_intergenic->length();

        #Write sequence to file. 

        if ($current_frame eq $last_gene_frame) {
            $counter[0] +=1;
            $outfile_intergenic->write_seq($output_intergenic);

        }elsif (($current_frame eq '-') and ($last_gene_frame eq '+')){
            $counter[1] +=1;
            $outfile_intergenic->write_seq($output_intergenic);

        }elsif (($current_frame eq '+') and ($last_gene_frame eq '-')){
            $counter[2] +=1;
            #do not write this to the file, since it contains no regulatory sequences.
        }else{
            die $current_frame . $last_gene_frame ."\n";
        }

        #Prepare next intergenic sequence, update counters.
        $last_gene_end = $gene_end;
        $last_gene_frame = $current_frame;
        $last_gene_id = $array[0];

    } elsif ($type eq 'gene') {
        #First gene.
        $last_gene_end = $array[4];
        $last_gene_frame = $array[6];
        $last_gene_id = $array[0];
        print("Found first gene.\n")
    }
}

print join(", ", @counter) . " total length: " . $total_length . "\n" ;
close GFF;
