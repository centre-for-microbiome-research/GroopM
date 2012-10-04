                                                                             
          .d8888b.                                    888b     d888          
         d88P  Y88b                                   8888b   d8888          
         888    888                                   88888b.d88888          
         888        888d888 .d88b.   .d88b.  88888b.  888Y88888P888          
         888  88888 888P"  d88""88b d88""88b 888 "88b 888 Y888P 888          
         888    888 888    888  888 888  888 888  888 888  Y8P  888          
         Y88b  d88P 888    Y88..88P Y88..88P 888 d88P 888   "   888          
          "Y8888P88 888     "Y88P"   "Y88P"  88888P"  888       888          
                                             888                             
                                             888                             
                                             888                             

Overview
=========

GroopM is a metagenomic binning toolset. It leverages spatio-temoral 
dynamics to accurately (and almost automatically) extract genomes 
from multi-sample metagenomic datasets.

GroopM is largely parameter-free. Use: groopm -h for more info.


Installation
=========

Should be as simple as

    pip install GroopM

Data preparation and running GroopM
=========

Before running GroopM you need to prep your data. A typical workflow looks like this:

    1. Produce NGS data for your environment across mutiple (3+) samples (spearated spatially or temporally or both).
    2. Co-assemble your reads using Velvet or similar.
    3. For each sample, map the reads against the co-assembly. GroopM needs sorted indexed bam files. If you have 3 samples then you will produce 3 bam files. I use BWA / Samtools for this.
    4. Take your co-assembled contigs and bam files and load them into GroopM using 'groopm parse' saveName contigs.fa bam1.bam bam2.bam...
    5. Keep following the GroopM workflow. See: groopm -h for more info.

Licence and referencing
=========

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/GroopM

This software is currently unpublished but a manuscript is being prepared. Please contact me at m_dot_imelfort_at_uq_dot_edu_dot_au for more information about referencing this software.

Copyright © 2012 Michael Imelfort. See LICENSE.txt for further details.
