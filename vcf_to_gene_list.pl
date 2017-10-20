#!/usr/bin/env perl
# 
# Author:       Mareike Herzog	
# Maintainer:   Mareike Herzog
# Created:      

use Carp;
use strict;
use warnings;
use Getopt::Long;

###############################
##							 ##
## CSQs we are interested in ## 
##							 ##
###############################
my @csqs = ('stop_gained', 'missense_variant','stop_lost','frameshift_variant','initiator_codon_variant','splice_donor_variant','splice_acceptor_variant','inframe_insertion','inframe_deletion','splice_region_variant');
my @samples;
my $aminoacid = 'NA'; my $nt_change = 'NA';

#############
##  INPUT  ##
#############
my ($input);

GetOptions
(
    'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];

my $ifh;
my $check = 0;
if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}

######################
##  Bulk of script  ##
######################
print qq[TYPE\tCHROM\tPOS\tGENE\tCSQ\tAminoAcid\tNtChange\tNumberofHOMS\tAllOtherCSQs\t];
while( my $l = <$ifh> )
{
    my $m = $l;
    $check = 0;
    chomp( $l );
    #skip all the header lines
    next if( $l =~ /^#/ && $l !~ /^#CHROM/);
    
    #split each line into its columns
    my @s = split( /\t/, $l );
    
    #get all the samples and print them in the header
    if( $l =~ /#CHROM/ ){for(my $i=9;$i<@s;$i++){$samples[$i]=$s[$i];print("$s[$i]\t");next;}print("\n");}
    
    #########################
    ## Find the right CSQs ##
    #########################
    #we go through each consequence
    foreach my $csq(@csqs)
    {
        
        #and check that the line matches this consequence
        if( $l =~ /$csq/i ) 
        {
    		
    		#if it does we split the info column into its components
            my @s1 = split( /;/, $s[ 7 ] );
            
            #and loop through each tag, we skip anything that does not start with CSQ
            foreach my $tag(@s1)
            {
                next unless $tag =~ /^(CSQ=)(.+)/;
                
                #the CSQ tag will be split on the + 
                my @s2 = split( /\+/, $2 );
                #and then we loop through the different consequences
                foreach my $csq(@s2)
                {
                    foreach my $csq_del(@csqs)
                    {
                        if( $csq =~ /$csq_del/i && $check ==0)
                        {
                            
                            #####################################
                            ## get number of samples affected  ##
                            #####################################
                            my $homs = 0; my $c = 0;
                            for(my $i=0;$i<@s;$i++){
                            my @genoty = split(/:/, $s[$i]);
                            if($genoty[0]=~/1\/1/){$homs++;}}
                            
                            #split the consequence on the :
                            my @s3 = split( /\:/, $csq );
                            
                            #################################
                            ##  Get the amino acid changes ##
                            #################################
                            $nt_change = "$s[3]".'>'."$s[4]";
                            if ( $csq =~ /missense/i )
                        	{
                        		#wanted: C911Y;
                        		#YLR442C[0]:YLR442C[1]:missense_variant[2]:2732[3]:911[4]:C>Y[5]
                                my @s5 = split( /\>/, $s3[5]);
                                $aminoacid = ' '.$s5[0].$s3[4].$s5[1].' ';
                                my $remainder = $s3[3] % 3;
                        		my $codon;
                        		if ($remainder == 0){$codon=3;}
                        		if ($remainder == 1){$codon=1;}
                        		if ($remainder == 2){$codon=2;}
                        		#make REF>ALT:Codon
                        		$nt_change = "$s[3]".'>'."$s[4]".':'."$codon";
                            }
                            if ( $csq =~ /initiator_codon_variant/i )
                        	{
                        		#wanted: M1I;
                        		#YIL129C[0]:YIL129C[1]:initiator_codon_variant[2]:3[3]:1[4]:M>I[5]
                        		my @s5 = split( /\>/, $s3[5]);
                        		$aminoacid = ' '.$s5[0].$s3[4].$s5[1].' ';
                        		my $remainder = $s3[3] % 3;
                        		my $codon;
                        		if ($remainder == 0){$codon=3;}
                        		if ($remainder == 1){$codon=1;}
                        		if ($remainder == 2){$codon=2;}	
                        		#make REF>ALT:Codon
                        		$nt_change = "$s[3]".'>'."$s[4]".':'."$codon";
                            }
                            if ( $csq =~ /frameshift_variant/i )
                        	{
                        		#wanted: range
                        		#YLL066W-B[0]:YLL066W-B[1]:frameshift_variant,feature_truncation[2]:59[3]:20[4]
                                $aminoacid = 'FS@'.$s3[4].' ';
                            }
                            if ( $csq =~ /stop_gained/i )
                            {
                                if ($m =~ /INDEL/) {
                                    #delta(number-1)
                                    #YLR442C[0]:YLR442C[1]:stop_gained[2]:2733-2775[3]:911-20[4]
                                    my @new_num = split ( /\-/, $s3[4] );
                                    my $number=$new_num[0]-1;
                                    $aminoacid = ' '.'Δ(FS)'.$number.' ';}
                                
                                else{
                                    #delta(number-1)
                                    #YLR442C[0]:YLR442C[1]:stop_gained[2]:2733[3]:911[4]
                                    my $number=$s3[4]-1;
                                    $aminoacid = ' '.'Δ'.$number.' ';}
                                    my $remainder = $s3[3] % 3;
                        			my $codon;
                        			if ($remainder == 0){$codon=3;}
                        			if ($remainder == 1){$codon=1;}
                        			if ($remainder == 2){$codon=2;}
                        			#make REF>ALT:Codon
                        			$nt_change = "$s[3]".'>'."$s[4]".':'."$codon";
                            }
                            if ( $csq =~ /synonymous/i )
                        	{
                        		#wanted £P91:C>T:2£
                        		#YDR420W[0]:YDR420W[1]:synonymous_variant[2]:1533[3]:511[4]:P>P[5]
                        		my @s5;
                        		if (defined $s3[5]){
                        		@s5 = split( />/, $s3[5]);}
                        		my $remainder = $s3[3] % 3;
                        		my $codon;
                        		if ($remainder == 0){$codon=3;}
                        		if ($remainder == 1){$codon=1;}
                        		if ($remainder == 2){$codon=2;}
                        		if (defined $s3[5]){
                        		$aminoacid = ' '.$s5[0].$s3[4].':'.$s[3].'>'.$s[4].':'.$codon.' ';}
                        		else{$aminoacid = ' '.$s3[4].':'.$s[3].'>'.$s[4].':'.$codon.' ';}
                        		#make REF>ALT:Codon
                        		$nt_change = "$s[3]".'>'."$s[4]".':'."$codon";
                            }
                           
                            #############
                            ##  PRINT  ##
                            #############
                            
                            #if( $homs > 0)
                            #{
                            if( $l=~/INDEL/){print qq[INDEL\t];}else{print qq[SNP\t];}
                            print qq[$s[0]\t$s[1]\t$s3[1]\t$s3[2]\t$aminoacid\t$nt_change\t$homs\t$tag\t];
                            for(my $i=9;$i<@s;$i++){$samples[$i]=$s[$i];print("$s[$i]\t");next;}print("\n");
                            $check = 1;
                            $aminoacid = 'NA'; $nt_change = 'NA';
                            
                            #}
                        }
                    }
                }
            }
        }
    }
}
close( $ifh );