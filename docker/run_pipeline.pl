#!/usr/bin/perl
# $Revision: 1.2 $
# $Date: 2021/05/26 $
# $Id: run_pipeline.pl $

# 
# Copyright 2014-2021, Michaël Bekaert <michael.bekaert@stir.ac.uk>
#
# This file is part of Vibrio-TNseq.
#
# Vibrio-TNseq is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Vibrio-TNseq is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License v3
# along with Vibrio-TNseq. If not, see <http://www.gnu.org/licenses/>.
#

use strict;
use warnings;
use Getopt::Long;
use Cwd qw/ abs_path /;

#----------------------------------------
our $VERSION = 1.2;

#----------------------------------------
my ($verbose, $threads, $outdir, $stats, $end, $gff, $cgview, $topology, $png, $svg, $database, $readfiles) = (0, 1, abs_path('.'), 0, 0, 0, 0, 'circular', 0, 0, '/databases/vibrio');
GetOptions(
           'i|infolder|readfiles:s' => \$readfiles,
           'stats:+'                => \$stats,
           'end:i'                  => \$end,
           'gff:s'                  => \$gff,
           'cgview:+'               => \$cgview,
           'topology:s'             => \$topology,
           'png:+'                  => \$png,
           'svg:+'                  => \$svg,
           'o|output:s'             => \$outdir,
           't|threads:+'            => \$threads,
           'v|verbose+'             => \$verbose
          );
if (defined $readfiles && -d $readfiles && defined $topology && $topology =~ m/(circular|linear)/sxm && $outdir && (-d abs_path($outdir) || !-e abs_path($outdir)) && defined $threads && int($threads) > 0)
{
    mkdir abs_path($outdir) if (!-e abs_path($outdir));
    opendir(my $dh, $readfiles) || die "Can't opendir $readfiles: $!";
    my @compressed = grep { /\.fastq(\.gz)?$/ && -f "$readfiles/$_" } readdir($dh);
    closedir $dh;
    foreach my $dataset (@compressed)
    {
        if ($dataset =~ m/^(.*)\.fastq(\.gz)?$/)
        {
            print {*STDERR} '> Processing ' . $dataset . "\n";
            print {*STDERR} '   * cutadapt -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTTCAGAGTTCTACAGTCCGACGATCACAC -a TAACAGGTTGGATGATAAGTCCCCGGTCTCTGTCTCTTATACACATCTCCGAGCCCACGAGAC -O 3 -m 10 -M 18 -e 0.15 --times 2 --trimmed-only -o ',
                   abs_path($outdir . q{/} . $1 . '.fastq'), q{ },
                   abs_path($readfiles . q{/} . $dataset), ' >',
                   abs_path($outdir . q{/} . $1 . '.log'), "\n" if ($verbose);
            system(  'cutadapt -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTTCAGAGTTCTACAGTCCGACGATCACAC -a TAACAGGTTGGATGATAAGTCCCCGGTCTCTGTCTCTTATACACATCTCCGAGCCCACGAGAC -O 3 -m 10 -M 18 -e 0.15 --times 2 --trimmed-only -o '
                   . abs_path($outdir . q{/} . $1 . '.fastq') . q{ }
                   . abs_path($readfiles . q{/} . $dataset) . ' >'
                   . abs_path($outdir . q{/} . $1 . '.log'));
            print {*STDERR} '   * bowtie2',
                   (int($threads) > 1 ? ' -p ' . int($threads) : q{}),
                   ' --no-1mm-upfront --end-to-end --very-fast -x ',
                   $database, ' -U ',
                   abs_path($outdir . q{/} . $1 . '.fastq'), ' -S ',
                   abs_path($outdir . q{/} . $1 . '.sam'), ' 2>>',
                   abs_path($outdir . q{/} . $1 . '.log'), "\n" if ($verbose);     
            system(  'bowtie2'
                   . (int($threads) > 1 ? ' -p ' . int($threads) : q{})
                   . ' --no-1mm-upfront --end-to-end --very-fast -x '
                   . $database . ' -U '
                   . abs_path($outdir . q{/} . $1 . '.fastq') . ' -S '
                   . abs_path($outdir . q{/} . $1 . '.sam') . ' 2>>'
                   . abs_path($outdir . q{/} . $1 . '.log'));
            system('rm ' . abs_path($outdir . q{/} . $1 . '.fastq'));
            print {*STDERR} '   * sam_to_map.pl -sam ',
                   abs_path($outdir . q{/} . $1 . '.sam'),
                   (defined $gff      ? ' --gff ' . $gff           : q{}),
                   (int($cgview) > 0  ? ' --cgview 1'              : q{}),
                   (defined $topology ? ' --topology ' . $topology : q{}),
                   (int($png) > 0     ? ' --png 1'                 : q{}),
                   (int($svg) > 0     ? ' --svg 1'                 : q{}),
                   (int($stats) > 0   ? ' --stats 1  --summary 1'  : q{}),
                   (int($end) > 0 && int($end) < 100 ? ' --end ' . $end : q{}),
                   (int($verbose) > 0 ? ' -v '                     : q{}) . ' -o ',
                   abs_path($outdir . q{/} . $1), ' >>',
                   abs_path($outdir . q{/} . $1 . '.log'), "\n" if ($verbose);
            system(  'sam_to_map.pl -sam '
                   . abs_path($outdir . q{/} . $1 . '.sam')
                   . (defined $gff      ? ' --gff ' . $gff           : q{})
                   . (int($cgview) > 0  ? ' --cgview 1'              : q{})
                   . (defined $topology ? ' --topology ' . $topology : q{})
                   . (int($png) > 0     ? ' --png 1'                 : q{})
                   . (int($svg) > 0     ? ' --svg 1'                 : q{})
                   . (int($stats) > 0   ? ' --stats 1 --summary 1'   : q{})
                   . (int($end) > 0 && int($end) < 100 ? ' --end ' . $end : q{})
                   . (int($verbose) > 0 ? ' -v '                     : q{}) . ' -o '
                   . abs_path($outdir . q{/} . $1) . ' >>'
                   . abs_path($outdir . q{/} . $1 . '.log'));
        }
    }
}
else
{
    print {*STDOUT}
      "Usage: $0 --infolder <compressed raw read file folder> [options]\n\n    --infolder   path to the compressed raw read file folder. [mandatory]\n    --gff        provide gene annotations and output a GFF file\n    --end        ignore the insert if is is in the last n%. [default 0/off]\n    --stats      provides coverage statistics. [default 0/off]\n    --cgview     output a cgview XML file. [default 0/off]\n    --topology   specify the chromosome topology circular or linear for gview. [default circular]\n    --png        run gview and generate an PNG image file. [default 0/off]\n    --svg        run gview and generate an SVG image file. [default 0/off]\n    --output     directory to use when generation outputs. [default .]\n    --verbose    become very chatty. [default 0/off]\n\n";
    exit 1;
}
