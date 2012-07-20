#!/usr/bin/perl
###############################################################################
#
#    cleanAssem.pl
#    
#    Sort out an assembly. Remove short contigs, sort by length and re-reference headers
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
use File::Spec;
use File::Basename;

#CPAN modules
use Bio::SeqIO;

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
# make the working dir
my $global_working_dir = File::Spec->rel2abs($global_options->{'working_dir'});
mkdir $global_working_dir;

# we need a sequence map for sorting strings by length
my $seq_map_ID = 1;
my %global_seq_map = ();
my %global_len_map = ();
my %global_cov_map = ();
my %contig_map = ();

# set the short sequence cut off
my $global_cut_off_len = overrideDefault(100, 'length');

if($global_options->{'outC'} eq $global_options->{'outA'} or $global_options->{'in'} eq $global_options->{'outA'} or $global_options->{'outC'} eq $global_options->{'in'})
{
    croak "**ERROR: Filenames must be unique\n";
}

# open our output files
my $global_cont_mapping_fh = openWrite($global_working_dir."/".$global_options->{'outC'});
my $global_clean_assem_fh = openWrite($global_working_dir."/".$global_options->{'outA'});

my $seqio = Bio::SeqIO->new(-file => $global_options->{'in'}, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) {
  # get the sequence
  my $string = $seq->seq;
  
  # reject short sequences
  next if(length $string < $global_cut_off_len);
  
  # store the string
  $global_len_map{$seq_map_ID} = length $string;
  $global_seq_map{$seq_map_ID} = fasta_cut($string);
  
  # get a list of headers
  my $header = $seq->id;
  print $global_cont_mapping_fh $global_options->{'tag'}."_$seq_map_ID\t$header\n";
  
  $seq_map_ID++;
}

# sort these guys by length
foreach my $con_id (sort {int($global_len_map{$b}) <=> int($global_len_map{$a})}  (keys %global_len_map))
{
    print $global_clean_assem_fh ">".$global_options->{'tag'}."_$con_id"."_$global_len_map{$con_id}\n";
    print $global_clean_assem_fh $global_seq_map{$con_id};
}

# clean up
close $global_cont_mapping_fh; 
close $global_clean_assem_fh;

######################################################################
# CUSTOM SUBS
######################################################################

sub fasta_cut {
    #-----
    # Cut up a fasta sequence 
    #
    my ($string) =  @_;
    $string = uc $string;
    $string =~ s/[^ACGT]/N/g;
    my $line_len = 80;
    my $return_str = "";
    my $len = length $string;
    my $start = 0;
    while($start < $len)
    {
        $return_str .= substr $string, $start, $line_len;
        $return_str .="\n";
        $start += $line_len;
    }
    return $return_str;
}

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    my @standard_options = ( "help|h+", "working_dir|w:s", "in|i:s", "outC|C:s", "outA|A:s", "length|l:i", "tag|t:s" );
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
    if(!exists $options{'in'} ) { printParamError ("You MUST supply an assembly to clean"); }
    if(!exists $options{'outA'} ) { printParamError ("Please specify where you would like the clean assembly saved"); }
    if(!exists $options{'outC'} ) { printParamError ("Please specify where you would like the contig references saved"); }
    if(!exists $options{'tag'} ) { printParamError ("You need to specify a tag"); }
    if(!exists $options{'working_dir'} ) { printParamError ("Please supply a working directory"); }    
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

    cleanAssem.pl

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

   Sort out an assembly. Remove short contigs, sort by length and re-reference headers

=head1 SYNOPSIS

    cleanAssem.pl -in|i FILE -outC|C FILE -outA|A FILE -tag|t STRING -working_dir|w DIR

      -in|i FILE                   Assembly file to clean
      -outC|C FILE                 File to print contig references
      -outA|A FILE                 File to print cleaned assembly
      -tag|t STRING                A UNIQUE tag for this cleaned assembly
      -working_dir -w DIR          A place to put all the files
      [-length|l INT]              Reject all contigs shorter than this length (default: 100bp)
      [-help -h]                   Displays basic usage information
         
=cut
