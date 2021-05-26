#!/usr/bin/perl
# $Revision: 1.2 $
# $Date: 2021/05/26 $
# $Id: sam_to_map.pl $

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
sub print_fit
{
    my @tab = @_;
    if (scalar @tab >= 20)
    {
        return ';T0=' . $tab[1] . ';T1ctl=' . $tab[2] . ';T2ctl=' . $tab[4] . ';T1exp=' . $tab[3] . ';T2exp=' . $tab[5] . ';T1_pvalue=' . $tab[11]  . ';T2_pvalue=' . $tab[12] . ';sig=' . $tab[13] . ';shape=' . $tab[14] . ';Gctl=' . $tab[16] . ';Gexp=' . $tab[18] . ';Fctl=' . $tab[19] . ';Fexp=' . $tab[20];
    }
    return q{};
}
my ($verbose, $threads, $outdir, $summary, $stats, $start, $end, $cgview, $topology, $png, $svg, $thickness, $gview, $gff, $sam, $fitness) = (0, 1, abs_path('./output'), 0, 0, 0, 0, 0, 'circular', 0, 0, 0.0002, '/usr/local/bin/gview.jar');
GetOptions('sam=s' => \$sam, 'start:i' => \$start, 'end:i' => \$end, 'summary:+' => \$summary, 'stats:+' => \$stats, 'gff:s' => \$gff, 'fitness:s' => \$fitness, 'cgview:+' => \$cgview, 'topology:s' => \$topology, 'png:+' => \$png, 'svg:+' => \$svg, 'o|output:s' => \$outdir, 't|threads:+' => \$threads, 'v|verbose+' => \$verbose);
if (defined $sam && -r $sam && defined $topology && $topology =~ m/(circular|linear)/sxm && $outdir && (-d abs_path($outdir) || !-e abs_path($outdir)) && defined $threads && int($threads) > 0)
{
    my (%loc, %genes, %chrs, %fit);
    my ($header, $distribution, $with, $without) = (q{}, $sam, 'Genes with TN', 'Genes without TN');
    mkdir abs_path($outdir) if (!-e abs_path($outdir));
    if (defined $gff && -r $gff && open(my $in, q{<}, $gff))
    {
        print {*STDERR} "> Accessing annotation...\n" if ($verbose);
        while (<$in>)
        {
            next if (m/^#/);
            chomp;
            my @data = split /\t/;
            if (scalar @data >= 9 && $data[2] eq 'gene' && int($data[3]) > 0 && int($data[4]) > int($data[3]))
            {
                my $name = $data[0];
                $name = $1 if ($data[0] =~ m/ref\|([^\|]+)\.\d+\|/);
                $data[1] = 0;
                $data[5] = [];
                @{$genes{$name}{int($data[3])}} = @data;
            }
            elsif (scalar @data >= 9 && $data[2] eq 'CDS' && int($data[3]) > 0 && int($data[4]) > 0 && $data[8] =~ m/product=([^;]+);/)
            {
                my $product = $1;
                my $name = $data[0];
                $name = $1 if ($data[0] =~ m/ref\|([^\|]+)\.\d+\|/);
                if (exists $genes{$name}{int($data[3])}) { $genes{$name}{int($data[3])}->[8] .= ';product=' . $product; }
            }
        }
        close $in;
    }
    if (defined $fitness && -r $fitness && open(my $in, q{<}, $fitness))
    {
        print {*STDERR} "> Accessing fitness statistics...\n" if ($verbose);
        while (<$in>)
        {
            next if (m/^#/);
            chomp;
            my @data = split /\t/;
            if (scalar @data >= 20 && $data[0] =~ m/^(.+)\_(\d+)$/) { @{$fit{$1}{$2}} = @data; }
        }
        close $in;
    }
    if (open my $in, q{<}, $sam)
    {
        print {*STDERR} "> Accessing alignments...\n" if ($verbose);
        while (<$in>)
        {
            chomp;
            my @data = split /\t/;
            if (scalar @data >= 12 && (int($data[1]) == 0 || int($data[1]) == 16) && int($data[4]) > 20)
            {
                # 0 => forward read (add length), 16, reverse read
                my $tn = int($data[3]);
                if (int($data[1]) == 0) { $tn += length($data[9]); }
                push @{$loc{$data[2]}{$tn}}, [@data];
            }
            elsif (scalar @data >= 3 && $data[0] eq '@SQ' && $data[2] =~ m/^LN:(\d+)/)
            {
                my $len = $1;
                if ($data[1] =~ m/^SN:(.*)$/) { $chrs{$1} = $len; }
            }
        }
        close $in;
    }
    if (%loc)
    {
        my ($count, $chrcount) = (0, 0);
        print {*STDERR} '> Calculating distritution on ', (scalar keys %loc), " chromosomes/plamids\n" if ($verbose);
        open my $sumstat, q{>}, abs_path($outdir . '/tn.coverage.tsv') if (int($stats) > 0 && int($summary) > 0);
        print {$sumstat} "TN\tcoverage\tdistance_next\n" if (defined $sumstat);
        foreach my $chr (sort keys %loc)
        {
            my $name = $chr;
            $name = $1 if ($chr =~ m/ref\|([^\|]+)\.\d+\|/);
            $header .= "\t" . $name;
            my $last = 0;
            $chrcount++;
            open my $outtn, q{>}, abs_path($outdir . '/tn' . $chrcount . '.gff');
            print {$outtn} "##gff-version 3\n";
            open my $outstat, q{>}, abs_path($outdir . '/chr' . $chrcount . '.coverage.tsv') if (int($stats) > 0);
            print {$outstat} "TN\tcoverage\tdistance_next\n" if (defined $outstat);
            open my $out, q{>}, abs_path($outdir . '/chr' . $chrcount . '.xml') if (int($cgview) > 0);
            print {$out} '<?xml version="1.0" encoding="ISO-8859-1"?>' . "\n"
              . '<cgview backboneRadius="300" backboneColor="rgb(102,102,102)" backboneThickness="2" featureSlotSpacing="4" labelLineLength="60" labelPlacementQuality="best" labelLineThickness="1" rulerPadding="14" tickThickness="1" arrowheadLength="4" rulerFont="SansSerif, plain, 8" rulerFontColor="rgb(0,0,0)" labelFont="SansSerif, plain, 10" isLinear="'
              . ($topology eq 'linear' ? 'true' : 'false')
              . '" minimumFeatureLength="0.2" globalLabel="', (int($cgview) > 1 ? 'true' : 'false'),
              '" moveInnerLabelsToOuter="false" featureThickness="8" tickLength="5" useInnerLabels="true" shortTickColor="rgb(0,51,0)" longTickColor="rgb(0,51,0)" zeroTickColor="rgb(0,51,0)" showBorder="true" borderColor="black" backgroundColor="white" tickDensity="0.5" sequenceLength="'
              . $chrs{$chr}
              . '" title="', $chr, '" height="1200" width="900"><featureSlot strand="direct" showShading="false">', "\n"
              if (defined $out);
            $distribution .= "\t" . (scalar keys %{$loc{$chr}});

            foreach my $pos (sort { $a <=> $b } keys %{$loc{$chr}})
            {
                $count++;
                my $flagin;
                if (%genes && exists $genes{$name})
                {
                    foreach my $gene (sort { $a <=> $b } keys %{$genes{$name}})
                    {
                        my ($spadding,$epadding) = (0, 0);
                        if (defined $end && int($end) > 0 && int($end) < 100) { $epadding = int(($genes{$name}{$gene}->[4] - $genes{$name}{$gene}->[3]) * ($end / 100)); }
                        if (defined $start && int($start) > 0 && int($start) < 100) { $spadding = int(($genes{$name}{$gene}->[4] - $genes{$name}{$gene}->[3]) * ($start / 100)); }

                        if ($genes{$name}{$gene}->[6] eq '+')
                        {
                            next if ($genes{$name}{$gene}->[4] - $epadding < $pos);
                            if ($genes{$name}{$gene}->[3] + $spadding <= $pos && $genes{$name}{$gene}->[4] >= $pos - $epadding)
                            {
                                $genes{$name}{$gene}->[1]++;
                                # push @{$genes{$name}{$gene}->[5]}, $fit{$name}{$pos} if(exists $fit{$name}{$pos});
                                $flagin = $genes{$name}{$gene}->[8];
                            }
                            # last;
                        }
                        # elsif ($genes{$name}{$gene}->[6] eq '-')
                        else
                        {
						    next if ($genes{$name}{$gene}->[4] - $spadding < $pos);
						    if ($genes{$name}{$gene}->[3] + $epadding <= $pos && $genes{$name}{$gene}->[4] >= $pos - $spadding)
                            {
                                $genes{$name}{$gene}->[1]++;
                                # push @{$genes{$name}{$gene}->[5]}, $fit{$name}{$pos} if(exists $fit{$name}{$pos});
                                $flagin = $genes{$name}{$gene}->[8];
                            }
                            # last;
                        }
                        last;
                    }
                }
                print {$outstat} $name, q{_}, $pos, "\t", (scalar @{$loc{$chr}{$pos}}), "\t", ($last > 0 ? int($pos - $last) : q{}), "\n" if (defined $outstat);
                print {$sumstat} $name, q{_}, $pos, "\t", (scalar @{$loc{$chr}{$pos}}), "\t", ($last > 0 ? int($pos - $last) : q{}), "\n" if (defined $sumstat);
                $last = $pos;
                print {$outtn} $name, "\t-\tmisc_difference\t", $pos, "\t", ($pos + 1), "\t.\t", (int($loc{$chr}{$pos}[0][1]) == 16 ? q{-} : q{+}), "\t.\tnote=Insertion Site;coverage=", (scalar @{$loc{$chr}{$pos}}),
                  (defined $flagin ? (defined $fitness ? q{;} . $flagin : q{}) . ';gene=true' : q{}), (exists $fit{$name}{$pos} ? print_fit(@{$fit{$name}{$pos}}) : q{}), "\n";
                if (defined $out)
                {
                    my $start = $pos - ($chrs{$chr} * $thickness);
                    my $end = $pos + 1 + ($chrs{$chr} * $thickness);
                    if ($start < 1)         { $start = $chrs{$chr} + $start; }
                    if ($end > $chrs{$chr}) { $end   = $end - $chrs{$chr} + 1; }
                    print {$out} '<feature color="', (defined $flagin ? 'rgb(0,160,224)' : 'rgb(163,218,246)'), '" decoration="arc" label="', $chr, " ", $pos, " ", (scalar @{$loc{$chr}{$pos}}), '" showLabel="', (int($cgview) > 1 ? 'true' : 'false'),
                      '"><featureRange start="', int($start), '" stop="', int($end), '" /></feature>', "\n";
                }
            }
            close $outtn;
            print {$out} '</featureSlot>' if (defined $out);
            if (%genes && exists $genes{$name})
            {
                my (@glimmer_direct, @glimmer_reverse);
                my ($countwith, $countwithout) = (0, 0);
                foreach my $gene (sort { $a <=> $b } keys %{$genes{$name}})
                {
                    if ($genes{$name}{$gene}->[3] < 1)         { $genes{$name}{$gene}->[3] = $chrs{$chr} + $genes{$name}{$gene}->[3]; }
                    if ($genes{$name}{$gene}->[4] > $chrs{$chr}) { $genes{$name}{$gene}->[4]   = $genes{$name}{$gene}->[4] - $chrs{$chr} + 1; }
                    
                    if ($genes{$name}{$gene}->[6] eq '+')
                    {
                        push @glimmer_direct,
                            '<feature color="'
                          . ($genes{$name}{$gene}->[1] > 0 ? 'rgb(0,149,68)' : 'rgb(169,211,176)')
                          . '" decoration="arc"><featureRange start="'
                          . $genes{$name}{$gene}->[3]
                          . '" stop="'
                          . $genes{$name}{$gene}->[4]
                          . '" /></feature>';
                    }
                    else
                    {
                        push @glimmer_reverse,
                            '<feature color="'
                          . ($genes{$name}{$gene}->[1] > 0 ? 'rgb(230,79,41)' : 'rgb(248,191,162)')
                          . '" decoration="arc"><featureRange start="'
                          . $genes{$name}{$gene}->[3]
                          . '" stop="'
                          . $genes{$name}{$gene}->[4]
                          . '" /></feature>';
                    }
                    if   ($genes{$name}{$gene}->[1] > 0) { $countwith++; }
                    else                                 { $countwithout++; }
                }
                print {$out} '<featureSlot strand="direct" showShading="false">',  "\n", join("\n", @glimmer_direct),  "\n", '</featureSlot>', "\n" if (defined $out && scalar @glimmer_direct);
                print {$out} '<featureSlot strand="reverse" showShading="false">', "\n", join("\n", @glimmer_reverse), "\n", '</featureSlot>', "\n" if (defined $out && scalar @glimmer_reverse);
                $with    .= "\t" . $countwith;
                $without .= "\t" . $countwithout;
            }
            print {$out} '<legend position="upper-center" font="SansSerif, plain, 20"><legendItem text="', $chr, '" textAlignment="center" /></legend>', "</cgview>\n" if (defined $out);
            close $out     if (defined $out);
            close $outstat if (defined $outstat);
            if (int($cgview) > 0 && int($png) > 0)
            {
                print {*STDERR} '   * ', 'java -Xmx4G -jar ', $gview, ' -l ', $topology, ' -W 2340 -H 1820 -f png -i ', abs_path($outdir . '/chr' . $chrcount . '.xml'), ' -o ', abs_path($outdir . '/chr' . $chrcount . '.png'), "\n" if ($verbose);
                if (system('java -Xmx4G -jar ' . $gview . ' -l ' . $topology . ' -W 2340 -H 1820 -f png -i ' . abs_path($outdir . '/chr' . $chrcount . '.xml') . ' -o ' . abs_path($outdir . '/chr' . $chrcount . '.png')) != 0)
                {
                    print {*STDERR} "WARNING: gview failed to execute: $!\n";
                }
            }
            if (int($cgview) > 0 && int($svg) > 0)
            {
                print {*STDERR} '   * ', 'java -Xmx4G -jar ', $gview, ' -l ', $topology, ' -W 2340 -H 1820 -f svg -i ', abs_path($outdir . '/chr' . $chrcount . '.xml'), ' -o ', abs_path($outdir . '/chr' . $chrcount . '.svg'), "\n" if ($verbose);
                if (system('java -Xmx4G -jar ' . $gview . ' -l ' . $topology . ' -W 2340 -H 1820 -f svg -i ' . abs_path($outdir . '/chr' . $chrcount . '.xml') . ' -o ' . abs_path($outdir . '/chr' . $chrcount . '.svg')) != 0)
                {
                    print {*STDERR} "WARNING: gview failed to execute: $!\n";
                }
            }
        }
        close $sumstat if (defined $sumstat);
        $header .= "\tTN";
        $distribution .= "\t" . $count;
        if (%genes && open my $out, q{>}, abs_path($outdir . '/genes.gff'))
        {
            print {*STDERR} "> Writing gene distribution...\n" if ($verbose);
            print {$out} "##gff-version 3\n";
            foreach my $chr (sort keys %genes)
            {
                foreach my $gene (sort { $a <=> $b } keys %{$genes{$chr}})
                {
                    print {$out} $genes{$chr}{$gene}->[0], "\t-\t", $genes{$chr}{$gene}->[2], "\t", $genes{$chr}{$gene}->[3], "\t", $genes{$chr}{$gene}->[4], "\t", $genes{$chr}{$gene}->[1], "\t", $genes{$chr}{$gene}->[6], "\t.\t",
                      $genes{$chr}{$gene}->[8], "\n";
                }
            }
            close $out;
        }
        print {*STDOUT} $header, "\n", $distribution, "\n", (%genes ? $with . "\n" . $without . "\n" : q{}), "\n";
    }
}
else
{
    print {*STDOUT}
      "Usage: $0 --sam <bowtie2 sam output> [options]\n\n    --sam        path to the SAM file to process must have SQ headers. [mandatory]\n    --gff        provide gene annotations and output a GFF3 file.\n    --fitness    provide fitness file to be incorporated [test only].\n    --start      ignore the insert if is is in the first n%. [default 0/off]    --end        ignore the insert if is is in the last n%. [default 0/off]\n    --stats      provides coverage statistics. [default 0/off]\n    --summary    summarised the stistics output. [default 0/off]\n    --cgview     output a cgview XML file. [default 0/off]\n    --topology   specify the chromosome topology circular or linear for gview. [default circular]\n    --png        run gview and generate an PNG image file. [default 0/off]\n    --svg        run gview and generate an SVG image file. [default 0/off]\n    --output     directory to use when generation outputs. [default ./output]\n    --verbose    become very chatty. [default 0/off]\n\n";
    exit 1;
}
