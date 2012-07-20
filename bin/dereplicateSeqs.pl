#!/usr/bin/env perl
###############################################################################
#
#    dereplicateSeqs.pl
#    
#    Combine the outputs of multiple runs of shatter.pl and dereplicate.
#
#    Copyright (C) Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use Carp;

#CPAN modules
#CPAN modules
use File::Spec;
use File::Basename;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Seq;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands 
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
printAtStart();
my $global_options = checkParams();

######################################################################
# CODE HERE
######################################################################
my $global_minimum_identity = overrideDefault(95,'global_minimum_identity');
my $global_length_cutoff = overrideDefault(95,'consumed_length');
$global_length_cutoff /= 100;
my $global_blast_threads = overrideDefault(1, 'blast_threads'); # the number of threads to give to blast

# make the working dir
my $global_working_dir = File::Spec->rel2abs($global_options->{'working_dir'});
mkdir $global_working_dir;

print "[$0] Preprocessing...\n";
# need a file for the combined shattered stuff
my $global_combined_file = $global_working_dir."/combined.fasta";

# first, combine all the files into one lump
my @shat_files = split(/,/, $global_options->{'infiles'});
my $file_list = join " ", @shat_files;
`cat $file_list > $global_combined_file`;
print "[$0] Dereplicating sequences in: $file_list\n";


print "[$0] Creating blast database...\n";
# now make a blast database
formatDB($global_combined_file);

print "[$0] Running blast, this make take some time...\n";
# now blast this sucker against itself
my $blast_file = blast($global_combined_file);

print "[$0] Identifying replicates [IDENT: $global_minimum_identity"."%, LEN: ".($global_length_cutoff*100)."%], this make take some time...\n";
# find all the replicated guys!
my %rep_list = findReps($blast_file);

print "[$0] Dereplicating input...\n";
# now run through the combined file and remove all the sequences which
# are replicates
my $out_file = $global_working_dir."/".$global_options->{'outfile'};
my $out_fh = openWrite($out_file);
my $total_seqs = 0;
my $kept_seqs = 0;
my $seqs = Bio::SeqIO->new(-file   => $global_combined_file,
                            -format => 'fasta'
                           );
while (my $seq = $seqs->next_seq) {
    $total_seqs++;
    if(!exists $rep_list{$seq->id})
    {
        $kept_seqs++;
        print $out_fh ">".$seq->id."\n";
        print $out_fh $seq->seq,"\n";
    }
}
close $out_fh;

# make the blast db for the dereplicated super set
print "[$0] Making dereplicated blast db\n";
formatDB($out_file);

# clean up
if(!exists $global_options->{'keep_temp_files'})
{
    removeFile($blast_file);
    removeFile($global_combined_file.".*");
    removeFile("formatdb.log");
}

print "[$0] DONE: Processed: $total_seqs and retained: $kept_seqs\n";

######################################################################
# CUSTOM SUBS
######################################################################
sub removeFile
{
    #-----
    # Wrap RM in our error code
    #
    my ($file) = @_;
    if(-e $file) {
        checkAndRunCommand("rm", {
                                 "" => $file
                                 }, WARN_ON_FAILURE);
    }
}

sub findReps
{
    #-----
    # Find all the replicated genes across the shattered files
    #
    my ($blast_file) = @_;
    my %gene_lengths = ();
    my %rep_list = ();
    my $in = new Bio::SearchIO(-format => 'blasttable', 
                               -file   => $blast_file);
                               
    # first, pick up some info we need
    while( my $result = $in->next_result ) 
    {
        my $q_name = $result->query_name;
        if(!exists $gene_lengths{$q_name})
        {
            my @header_fields = split(/_/, $q_name);
             $gene_lengths{$q_name} = (int($header_fields[$#header_fields - 1]) - int($header_fields[$#header_fields - 2]) + 1);
             
        }
        my $q_length = $gene_lengths{$q_name};
        while( my $hit = $result->next_hit ) 
        {
            my $h_name = $hit->name;
            next if($q_name eq $h_name);
            if(!exists $gene_lengths{$h_name})
            {
                my @header_fields = split(/_/, $h_name);
                 $gene_lengths{$h_name} = (int($header_fields[$#header_fields - 1]) - int($header_fields[$#header_fields - 2]) + 1);
                 
            }
            my $h_length = $gene_lengths{$h_name};
            while( my $hsp = $hit->next_hsp ) 
            {
                # now test if it meets our exacting requirements :)
                if ( $hsp->percent_identity >= $global_minimum_identity ) {
                    my $olap_length = $hsp->length('total'); 
                    if( ($olap_length >= $h_length) or 
                        ($olap_length >= $q_length) or 
                        ($olap_length/$h_length >= $global_length_cutoff) or 
                        ($olap_length/$q_length >= $global_length_cutoff) )
                    {
                        # replicate!
                        if($h_length < $q_length) 
                        {
                            $rep_list{$h_name} = 1;
                        }
                        elsif($q_length < $h_length) 
                        {
                            $rep_list{$q_name} = 1;
                        }
                        else 
                        {
                            # they may be the same length!
                            if($h_name lt $q_name)
                            {
                                $rep_list{$h_name} = 1;
                            }
                            else
                            {
                                $rep_list{$q_name} = 1;
                            }
                        } 
                    }
                }
            }  
        }
    }    
    return %rep_list;
}

sub blast
{
    #-----
    # Wrap the bare blast command
    # 
    my ($blast_db_name) = @_;
    my $blast_out = $blast_db_name.".blast";
    
    if(exists $global_options->{'legacy_blast'})
    {
        checkAndRunCommand("blastall", {
                                       -p => "blastn",
                                       -d => $blast_db_name,
                                       -i => $global_combined_file,
                                       -m => "8",
                                       -o => $blast_out
                                       }, DIE_ON_FAILURE);        
    }
    else
    {
        checkAndRunCommand("blastn", {
                                     -db => $blast_db_name,
                                     -query => $global_combined_file,
                                     -outfmt => "6",
                                     -out => $blast_out,
                                     -num_threads => $global_blast_threads
                                     }, DIE_ON_FAILURE);        
    }
    return $blast_out;
}

sub formatDB
{
    #-----
    # wrap the bare command - perhaps this should be handled by bioperl?
    #
    my ($fasta_file) = @_;
    if(exists $global_options->{'legacy_blast'})
    {
        checkAndRunCommand("formatdb", {
                                       -p => "F",
                                       -o => "T",
                                       -i => $fasta_file
                                       }, DIE_ON_FAILURE);        
    }
    else
    {
        checkAndRunCommand("makeblastdb", {
                                          -in => $fasta_file,
                                          -dbtype => "nucl",
                                          -parse_seqids => ""  
                                          }, DIE_ON_FAILURE);
    }
}

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ("help|h+", "infiles|f:s", "outfile|o:s", "working_dir|w:s", "blast_threads|b:i", "legacy_blast|L+", "keep_temp_files|k+", "minimum_identity|i:i", "consumed_length|l:i");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    if(!exists $options{'infiles'} ) { printParamError ("Need a set of input files to continue"); }
    if(!exists $options{'outfile'} ) { printParamError ("Need an output file name to continue"); }
    if(!exists $options{'working_dir'} ) { printParamError ("Please supply a working directory"); }
    if(exists $options{'blast_threads'} && exists $options{'legacy_blast'} ) { printParamError ("Multithreading only works with newer versions of Blast"); }
    #if(!exists $options{''} ) { printParamError (""); }

    return \%options;
}
sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #  
    my ($error) = @_;  
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
#                          -a => ""
#                          }, 
#                          WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    # 
    my ($cmd, $params, $failure_type) = @_;
    
    isCommandInPath($cmd, $failure_type);
    
    # join the parameters to the command
    my $param_str = join " ", map { $_ . " " . $params->{$_}} keys %{$params};
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    # make sure that all went well    
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : "  . $! . "\n";
        }
    }
}

######################################################################
# MISC

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    dereplicateSeqs.pl

=head1 COPYRIGHT

   copyright (C) Michael Imelfort

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

   This script takes a bunch of fasta files as input and generates a single file
   which contains a dereplicated subset of the input. It uses Blast to do the 
   matching, by default it looks for 95%+ overlap and 95%+ identity. Thus it has
   been developed with shorter sequences (Genes?) in mind.

    *** requires:
    
    Blast

=head1 SYNOPSIS

    dereplicateSeqs.pl -infiles|f FILE,FILE[,FILE,..] -working_dir|w DIR -outfile|o FILE
     
      -infiles -f FILE,FILE[,FILE,..]           The files we wish to combine and dereplicate
      -outfile -o FILE                          Write to this file
      -working_dir -w WORKING_DIR               A place to put all the files

   Parameters used for determining replication
   
      [-minimum_identity -i INT]                The minimum identity across this length [default: 95]      
      [-consumed_length -l INT]                 The amount of the smaller gene which must be present in the larger [default: 95]      
   
   Parameters used when running Blast   
      
      [-blast_threads -b NUM_THREADS]           The number of threads ot give to blast [default: 1]
      [-legacy_blast -L]                        Check this to use blast 2.2.22 etc... [default: expects 2.2.25]
                                                   NOTE: Incompatible with multiple threads
    Misc

      [-keep_temp_files -k]                     Keep the temp files? [default: delete files at the end]
      [-help -h]                                Displays basic usage information

=cut

