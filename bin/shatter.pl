#!/usr/bin/perl
###############################################################################
#
#    shatter.pl
#    
#    Run Glimmer on a set of cleaned contigs and create a new set of 
#    indexed orfs which are ready for fingerprinting...
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
use Cwd;
 
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
# check that the file exists
checkFileExists($global_options->{'in'});

# make the working dir
my $global_working_dir = File::Spec->rel2abs($global_options->{'working_dir'});
mkdir $global_working_dir;

if(exists $global_options->{'glimmer'})
{
    # run Glimmer3 (stupid g3 from scratch script!)
    isCommandInPath("g3-from-scratch.csh", DIE_ON_FAILURE);
    my $this_dir = cwd;
    cd $global_working_dir;
    `g3-from-scratch.csh $global_options->{'in'} $global_options->{'tag'}`;
    cd $this_dir;
    
    # make the glimmer output into gff3
    checkAndRunCommand("glim2gff3.pl", {
                                       -in => "$global_options->{'tag'}\.predict",
                                       -out => $global_working_dir."/$global_options->{'tag'}\.gff"
                                       }, DIE_ON_FAILURE);
    
    # remove all the glimmer files
    removeFile($global_options->{'tag'}.".detail");
    removeFile($global_options->{'tag'}.".icm");
    removeFile($global_options->{'tag'}.".longorfs");
    removeFile($global_options->{'tag'}.".predict");
    removeFile($global_options->{'tag'}.".train");
    
}
else
{
    # run prodigal
    checkAndRunCommand("prodigal", {
                                   -p => "meta",
                                   -i => $global_options->{'in'},
                                   -f => "gff",
                                   -o => $global_working_dir."/$global_options->{'tag'}\.gff"
                                   }, DIE_ON_FAILURE);
}

# make the gff3 file into a multiple fasta

checkAndRunCommand("gff2fasta.pl", {
                                   -gff => "$global_options->{'tag'}\.gff",
                                   -fasta => $global_options->{'in'},
                                   -out => $global_working_dir."/$global_options->{'tag'}\.shattered\.fasta",
                                   -include_nulls => "",
                                   -w => "0",
                                   -l => "150"
                                   }, DIE_ON_FAILURE);

# dun!

######################################################################
# CUSTOM SUBS
######################################################################
sub removeFile
{
    #-----
    # Wrap RM in our error code
    #
    my ($file) = @_;
    checkAndRunCommand("rm", {
                             "" => $file
                             }, WARN_ON_FAILURE);
}

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    my @standard_options = ( "help|h+", "in|i:s", "working_dir|w:s", "tag|t:s", "glimmer|g+" );
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
    if(!exists $options{'in'} ) { printParamError ("You MUST supply a file of contigs to parse"); }
    if(!exists $options{'tag'} ) { printParamError ("Please supply a tag"); }
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
    
    # check the command is in the path
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

    shatter.pl

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

   Run Prodigal or Glimmer on a set of cleaned contigs and create a new set of 
   indexed orfs which are ready for further analysis

    *** requires:
    
    Glimmer3 or Prodigal

=head1 SYNOPSIS

    shatter.pl -in|i CONTIGS -tag|t TAG -working_dir|w DIR 

      -in -i CONTIGS               Contig file to parse
      -tag -t TAG                  TAG passed throug to the glimmer program
      -working_dir -w WORKING_DIR  A place to put all the files

      [-glimmer -g]                Use glimmer
      [-help -h]                   Displays basic usage information
         
=cut
