#/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

if ($#ARGV != 2 ) {
	print "Usage perl struct_assign.pl <current working directory path> genelist.txt <path for SCRATCH1D run_SCRATCH-1D_predictors.sh>\n";
	exit;
}
##############################################################################


my($dir,$gene_list,$scratch_path) = @ARGV;

print "The arguments passed are :\n";
foreach my $temp(@ARGV) #prints arguments
{
print"$temp\n";
}
my @loci = ();

open(my $fh2,"<","$dir/$gene_list") or die "Can't open $gene_list\n"; #creates loci array from the list of loci provided
{
while(my $line = <$fh2>)
{
chomp($line);
push(@loci,$line);
}
}

print" The list of queries contain following names:\n"; #prints names of the queries in the genelist file 
foreach my $locus(@loci)
{
print"$locus\n";
}

print"You specified $dir as your current working directory in the arguments\n";

my @new =();
my @SS_class = ("SS_HELIX","SS_SHEET","SS_COIL");
my @ACC_class = ("EXPOSED","BURIED");
my $SS_seq =();
my $ACC_seq = ();


foreach my $locus (@loci)  # Get a file with all loci and add them here in for loop(it initializes names for all other files
{
my $file = $locus . ".nex";
my $file_fasta = "$locus.fasta";
my $updated_nex = $locus . ".updated" . ".txt";
my $new_file = $locus . "NEWMOD" . ".nex";
my $SS_file = $file_fasta . ".out" . ".ss";
my $ACC_file = $file_fasta . ".out" . ".acc";
###############Parsing only sequences from nexus file into a scratch file###############
open (my $FH,"$dir/$file") or die "Couldn't read file $!";
{
	while (my $line = <$FH>)
	{
		chomp($line);
		if ($line=~/Nexus/i || $line=~/Matrix/i || $line=~/\;/ || $line=~/\,/ || $line=~/format/i || $line eq /\s/) 
			{
			next;
			}
		else
			{
				push (@new,$line);
				open(my $FH1,">>","$dir/$updated_nex")or die "Couldn't read file $!";
				{
					print $FH1 "$line\n";
				}close($FH1);
			}
	}
}close($FH);

###########Create Consensus sequence############
print"Creating consensus sequence for $updated_nex\n";
create_consensus($updated_nex,$file_fasta,$locus);

############ Create a copy of original nexus file#####
print"Creating copy of original nexus file as $new_file\n";
system("cp $dir/$file $dir/$new_file");

############Adding assumption block in $new_file#########
open(my $FH,">>","$dir/$new_file") or die "Couldn't open $!\n";
{
print $FH "\nBEGIN ASSUMPTIONS;\n";
print "$new_file created \n";
}

#################Calling scratch1D on consensus sequence####
print "Calling Scratch1D on $file_fasta \n";
call_scratch($file_fasta);

#####################Adding Secondary structure result in charset block########################################
open(my $fh,"<","$dir/$SS_file");
{
while(my $line1 = <$fh>)
	{
		chomp($line1);
		if ($line1 =~ m/^>/)
		{
		next;
		}
		else
		{
		$SS_seq = $line1;
		}
	}
}
my @SS_array = split(//,$SS_seq);
open(my $FH2,">>","$dir/$new_file") or die "Couldn't open $!\n";
{
foreach my $class(@SS_class)
{
	print $FH2 "CHARSET $class = ";
	for(my $i=0;$i<=$#SS_array;$i++)
	{
		my $j = $i+1;
		if(($SS_array[$i] eq "H") and ($class =~ m/\bSS_HELIX\b/))
		{
		print $FH2 "$j ";
		}
		elsif(($SS_array[$i] eq "E") and ($class =~ m/\bSS_SHEET\b/))
		{
		print $FH2 "$j ";
		}
		elsif(($SS_array[$i] eq "C") and ($class =~ m/\bSS_COIL\b/))
		{
		print $FH2 "$j ";
		}
	}
print $FH2 ";\n";
}
###############################Adding Solvent accessibility result in Charset block#############
open(my $FH,"<","$dir/$ACC_file");
{
while(my $line1 = <$FH>)
	{
		chomp($line1);
		if ($line1 =~ m/^>/)
		{
		next;
		}
		else
		{
		$ACC_seq = $line1;
		}
	}
}

my @ACC_array = split(//,$ACC_seq);
foreach my $class1(@ACC_class)
{
	print $FH2 "CHARSET $class1 = ";
	for(my $i=0;$i<=$#ACC_array;$i++)
	{
		my $j = $i+1;
		if(($ACC_array[$i] eq "e") and ($class1 =~ m/\bEXPOSED\b/))
		{
		print $FH2 "$j ";
		}
		elsif(($ACC_array[$i] eq "-") and ($class1 =~ m/\bBURIED\b/))
		{
		print  $FH2 "$j ";
		}
	}
print $FH2 ";\n";
}

print $FH2 "END;\n";
}
}


##################################################Subroutines##################################################
sub create_consensus
{
#This subroutine creates a consensus sequence using MSA in the input file. It takes 3 arguments < sequences without nexus header from nexus file stored as <$updated_nex>,Output file name for consensus sequence<$file_fasta>,Gene name as for fasta header<$locus>
my $input_file = $_[0];
my $out_file = $_[1];
my $header = $_[2];
my %seq_hash =();
my $length =0;
my %w_hash =();
my %sum_weight =();
my %aa_weight =();

open(my $FH,"<","$dir/$input_file") or die "Couldn't read file $!\n";
{
	while(my $line = <$FH>)             #It creates a hash of array for the sequence file
		{
		chomp($line);
		#print"$line\n";
		my ($taxon,$seq) = split(/\s+/,$line);
		#print "$taxon\n";
		chomp($seq);
		my @sequences = split(//,$seq);
		for(my $j=0;$j<=$#sequences;$j++)
			{
			$seq_hash{$taxon}[$j] = $sequences[$j];
			}
			$length = scalar( @{ $seq_hash{$taxon} } );
			#print"The length of the array for $taxon is $length\n"; 
		}
}		
for(my $i=0;$i<=$length;$i++)               ##It counts different types of amino acid at a position 
{
my %aa_count =();
for my $temp(keys %seq_hash)
	{
		my $char1 = $seq_hash{$temp}[$i];
		chomp($char1);
		$aa_count{$char1}++;                  #creates a hash with aa and its count 
	}
				#print"position $i\n";
				my $diff_res = scalar(keys %aa_count);         #counts number of entries in array which is equal to different types of residue
				#print" The number of different residues at $i = $diff_res\n";
				for my $key (keys %seq_hash)                    #create a hash with position weights (1/r*s)
				{
					my $seq_char = $seq_hash{$key}[$i];
					if (exists $aa_count{$seq_char})
					{
					my $val = $aa_count{$seq_char};          #The number of times this residue occurred = no of sequences it is present
					my $weight = 1/($diff_res * $val);       # Weight of that residue as described by Henikoff 
					#print "$weight  for $seq_char position $i\n";
					$w_hash{$key}[$i] = $weight;              # create a matrix (hash of arrays) of weights for all residues in the sequence
					}
				}
}
foreach (keys %w_hash)            # get the sum of all individual aa weights. It gives the sequence weight
{
my $sum_seq = sum(@{$w_hash{$_}});
$sum_weight{$_} = $sum_seq;
#print "The sum for $_ is $sum_seq\n";
}

foreach (keys %seq_hash)            # Assign the sequence weight to individual aa in the sequence
{
 my $seq_w = $sum_weight{$_};
 #print" The weight for $_ is $seq_w\n";
 for(my $i=0;$i<=$length;$i++)
 	{
 	if( $seq_hash{$_}[$i] ne "-")
 		{
 		$aa_weight{$_}[$i] = $seq_w;        # A new hash stores sequence weight for each residue for overall sequence
 		}
 	else
 		{
 		$aa_weight{$_}[$i] = 0;
 		}
 	}
}
my @consensus = ();
for(my $i=0;$i<=$length;$i++)           # Find the maximum weight aa at the given position and identify that residue to create a consensus sequence
{
	my %aa_wc =();
	foreach (keys %seq_hash)
	{
		my $pos_char = $seq_hash{$_}[$i];
		my $pos_val = $aa_weight{$_}[$i];
		$aa_wc{$pos_char} = $aa_wc{$pos_char} + $pos_val;
	}
	foreach my $key (sort { $aa_wc{$b} <=> $aa_wc{$a} } keys %aa_wc)
	{
        #printf "%4d %s\n", $aa_wc{$key}, $key;
        push(@consensus,$key);
        last;
    }
	#print Dumper(\%aa_wc);
}
my $final_seq = join('',@consensus);

open(my $fh1,">>","$dir/$out_file") or die "Couldn't write to $!\n";
{
print $fh1 ">$header\n$final_seq";
}

#print"$final_seq\n";
#print Dumper(\%aa_weight);
#print Dumper(\%sum_weight);
}


sub call_scratch
{
my $input = shift;
my $out = $input . ".out";

print"input file is $input and calling scratch1D \n";
system("$scratch_path/run_SCRATCH-1D_predictors.sh $input $out 4");

return $out;
}