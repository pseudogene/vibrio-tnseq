#! /usr/bin/env python
# coding: utf8
#
# Copyright 2021 Michaël Bekaert. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#    2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
# THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
#
# This is an beta version of a python implementation of EL-ARTIST, adapted from MatLab
# Pritchard, J.R., Chao, M.C., Abel, S., Davis, B.M., Baranowski, C., Zhang, Y.J.,
# Rubin, E.J., Waldor, M.K., 2014. ARTIST: High-Resolution Genome-Wide Assessment of
# Fitness Using Transposon-Insertion Sequencing. PLoS Genet. 10, e1004782. 
# https://doi.org/10.1371/journal.pgen.1004782
#
#
import argparse
import math
import os
import re

import numpy as np
from Bio import SeqIO
from Bio import Seq
from hmmlearn import hmm
from sklearn.neighbors import KernelDensity


def parse_genomes(gff, fasta, feature='CDS', site='TA', verbose=False):
    print('>Accessing Fasta file...')
    chrlist = {}
    tnsites = {}
    with open(fasta, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            chrlist[record.id] = len(record.seq)
            tnsites[record.id] = [m.start() for m in re.finditer(site, str(record.seq))]
    if verbose:
        print(chrlist)

    print('>Accessing GFF/GTF file...')
    geneid = re.compile(r'gene_id "([^"]+)"')  # Assuming GFF2 / GTF
    nameid = re.compile(r'name "([^"]+)"')
    productid = re.compile(r'product "([^"]+)"')
    gene_biotype = re.compile(r'gene "([^"]+)"')
    tmp = {}
    with open(gff, 'r') as handle:
        line = handle.readline()
        if 'gff-version 3' in line:
            geneid = re.compile(r'ID=([^;]+)')
            nameid = re.compile(r'gene=([^;]+)')
            productid = re.compile(r'product=([^;]+)')
            gene_biotype = re.compile(r'gene_biotype=([^;]+)')
        # elif 'gff-version 2' in line:
        #     geneid = re.compile(r'gene_id "([^"]+)"')
        # else:
        #     print('Unknown format')
        #     exit(2)

        # Note: if a feature overlaps another one, it truncates the second feature.
        last = 0
        for line in handle:
            tabs = line.strip().split('\t')
            if len(tabs) >= 9 and tabs[2] == feature and tabs[0] in chrlist:
                m = geneid.search(tabs[8])
                if m is not None:
                    if tabs[0] not in tmp:
                        tmp[tabs[0]] = {}
                    if last + 1 < int(tabs[4]):
                        tmp[tabs[0]][max(int(tabs[3]), last + 1)] = {'start': max(int(tabs[3]), last + 1), 'end': int(tabs[4]), 'name': m.group(1), 'chromosome': tabs[0], 'strand': tabs[6], 'reads': 0, 'sites': 0, 'score': 0}
                        m = nameid.search(tabs[8])
                        if m is not None:
                            tmp[tabs[0]][max(int(tabs[3]), last + 1)]['gene'] = m.group(1)
                        m = productid.search(tabs[8])
                        if m is not None:
                            tmp[tabs[0]][max(int(tabs[3]), last + 1)]['product'] = m.group(1)
                        m = gene_biotype.search(tabs[8])
                        if m is not None:
                            tmp[tabs[0]][max(int(tabs[3]), last + 1)]['biotype'] = m.group(1)

                    last = int(tabs[4])

    print('>Order the features...')
    featurelist = {}
    for chromosome in sorted(tmp.keys()):
        last = None
        featurelist[chromosome] = {}
        for entry in sorted(tmp[chromosome].keys()):
            if last is None:
                featurelist[chromosome][1] = {'start': 1, 'end': tmp[chromosome][entry]['start'] - 1, 'name': 'IG_' + tmp[chromosome][entry]['name'], 'reads': 0, 'score': 0, 'sites': 0, 'chromosome': chromosome}
            else:
                featurelist[chromosome][tmp[chromosome][last]['end'] + 1] = {'start': tmp[chromosome][last]['end'] + 1, 'end': tmp[chromosome][entry]['start'] - 1, 'name': 'IG_' + tmp[chromosome][entry]['name'], 'reads': 0, 'score': 0, 'sites': 0, 'chromosome': chromosome}
            featurelist[chromosome][entry] = tmp[chromosome][entry]
            last = entry
        if chromosome in chrlist:
            featurelist[chromosome][tmp[chromosome][last]['end'] + 1] = {'start': tmp[chromosome][last]['end'] + 1, 'end': chrlist[chromosome], 'name': 'IG_' + chromosome + '_END', 'reads': 0, 'score': 0, 'sites': 0, 'chromosome': chromosome}
    del tmp

    print('>Map TN insertions site...')
    tnlist = {}
    for chromosome in sorted(tnsites.keys()):
        last = 1
        current = 1
        listfeatures = sorted(featurelist[chromosome].keys())
        tnlist[chromosome] = {}
        for entry in sorted(tnsites[chromosome]):
            entry += 1  # Fix coordinates issue between string and sequence
            while entry > current and len(listfeatures) > 0:
                last = current
                current = listfeatures.pop(0)
            if len(listfeatures) == 0 and entry > current:
                tnlist[chromosome][entry] = {'position': entry, 'chromosome': chromosome, 'feature_name': featurelist[chromosome][current]['name'], 'name_id': current}
                featurelist[chromosome][current]['sites'] += 1
            else:
                tnlist[chromosome][entry] = {'position': entry, 'chromosome': chromosome, 'feature_name': featurelist[chromosome][last]['name'], 'name_id': last}
                featurelist[chromosome][last]['sites'] += 1
    del tnsites
    return [featurelist, tnlist, chrlist]


def parse_reports(tn_report, tnlist):
    print('>Access TN insertions reports...')
    geneid = re.compile(r'coverage "(\d+)"')  # Assuming GFF2 / GTF
    tmp = {}
    with open(tn_report, 'r') as handle:
        line = handle.readline()
        if 'gff-version 3' in line:
            geneid = re.compile(r'coverage=(\d+);')
        # elif 'gff-version 2' in line:
        #     geneid = re.compile(r'coverage "(\d+)"')
        # else:
        #     print('Unknown format')
        #     exit(2)

        for line in handle:
            tabs = line.strip().split('\t')
            if len(tabs) >= 9 and 'Insertion Site' in tabs[8]:
                m = geneid.search(tabs[8])
                if m is not None:
                    if tabs[6] == '-' and int(tabs[3]) - 2 in tnlist[tabs[0]]:
                        tabs[3] = int(tabs[3]) - 2
                    if int(m.group(1)) >= 1 or int(tabs[3]) in tnlist[tabs[0]]:
                        if tabs[0] not in tmp:
                            tmp[tabs[0]] = {}
                        if int(tabs[3]) in tmp[tabs[0]]:
                            tmp[tabs[0]][int(tabs[3])]['count'] += int(m.group(1))
                        else:
                            tmp[tabs[0]][int(tabs[3])] = {'position': int(tabs[3]), 'count': int(m.group(1)), 'chromosome': tabs[0]}
    return tmp


def compact_reports(tn_window, chrlist, tnlist, featurelist):
    print('>Binning TN map...')
    binned = {}
    for chromosome in sorted(tnlist.keys()):
        last = 1
        current = 1
        listfeatures = sorted(featurelist[chromosome].keys())
        numwindows = math.floor(chrlist[chromosome] / tn_window) + 1
        binned[chromosome] = {}
        for bins in range(0, numwindows):
            entry = int(bins * tn_window + tn_window / 2)
            while entry > current and len(listfeatures) > 0:
                last = current
                current = listfeatures.pop(0)
            if len(listfeatures) == 0 and entry > current:
                binned[chromosome][bins] = {'position': entry, 'bin': bins, 'chromosome': chromosome, 'feature_name': featurelist[chromosome][current]['name'], 'name_id': current, 'count': 0, 'tn': []}
            else:
                binned[chromosome][bins] = {'position': entry, 'bin': bins, 'chromosome': chromosome, 'feature_name': featurelist[chromosome][last]['name'], 'name_id': last, 'count': 0, 'tn': []}
        for tn in tnlist[chromosome].keys():
            bins = math.floor(tn / tn_window)
            binned[chromosome][bins]['count'] += tnlist[chromosome][tn]['count']
            binned[chromosome][bins]['tn'].append(tn)

    return binned


def update_features(tnlist, featurelist):
    print('>Updating gene statistics...')

    for chromosome in sorted(tnlist.keys()):
        current = 1
        last = 1
        listfeatures = sorted(featurelist[chromosome].keys())
        for entry in sorted(tnlist[chromosome].keys()):
            while entry > current and len(listfeatures) > 0:
                last = current
                current = listfeatures.pop(0)
            if len(listfeatures) == 0 and entry > current:
                featurelist[chromosome][current]['score'] += 1
                featurelist[chromosome][current]['reads'] += tnlist[chromosome][entry]['count']
            else:
                featurelist[chromosome][last]['score'] += 1
                featurelist[chromosome][last]['reads'] += tnlist[chromosome][entry]['count']

    return featurelist


def checkseqdepth(tnlist, prefix, verbose=False):
    print('>Calculate coverage...')
    if verbose:
        print('chromosome\tnb_sites\tnb_hits\tfraction\treads')
    with open(prefix + '.Rtab', 'w') as scatter:
        scatter.write('chromosome\treads\tnb_hit\n')
    for chromosome in sorted(tnlist.keys()):
        total = 0
        hit = 0
        for bins in tnlist[chromosome].keys():
            total += tnlist[chromosome][bins]['count']
            if tnlist[chromosome][bins]['count'] > 0:
                hit += 1
        if total > 0:
            taprop = []
            for bins in tnlist[chromosome].keys():
                taprop.append((tnlist[chromosome][bins]['count']) / total)

            sensitivityanalysis = []
            t = math.floor(total / 300)
            for i in range(0, 300):
                x = np.random.multinomial(n=i * t + 1, pvals=taprop, size=10)
                x = ((x.T >= 1) * 1).sum(axis=0).mean()
                sensitivityanalysis.append([i * t + 1, x])

            with open(prefix + '.Rtab', 'a') as scatter:
                for i in sensitivityanalysis:
                    scatter.write(chromosome + '\t' + str(i[0]) + '\t' + str(i[1]) + '\n')

            if verbose:
                fractionhit = hit / len(tnlist[chromosome].keys())
                print(chromosome, len(tnlist[chromosome].keys()), hit, fractionhit, total, sep='\t')


def slidingwindow(tnwindows, featurelist, similation=1000, windows=5, pvalue=0.005):
    print('>Localise essential regions...')

    essentialregions = {}
    for chromosome in sorted(tnwindows.keys()):
        boots = np.random.uniform(high=len(tnwindows[chromosome]), size=[windows, similation])
        tmp_a = []
        essentialregions[chromosome] = {}

        for i in range(0, similation):
            total = 0
            numberhit = 0
            for j in range(0, windows):
                ind = int(boots[j, i])
                item = list(tnwindows[chromosome].values())[ind]
                total += item['count']
                if item['count'] > 0:
                    numberhit += 1
            tmp_a.append(total)
        bootstats = np.array(tmp_a)[:, np.newaxis]

        kde = KernelDensity(kernel='gaussian').fit(bootstats)
        for item, value in tnwindows[chromosome].items():
            sumwindowval = 0
            for j in range(0, windows):
                key = item + j
                if key >= len(tnwindows[chromosome]):  # Bacteria have circular chromosomes!
                    key -= len(tnwindows[chromosome])
                sumwindowval += tnwindows[chromosome][key]['count']
            if sumwindowval < 1:
                sumwindowval = 1
            values = np.asarray([value for value in range(0, int(sumwindowval))])
            values = values.reshape((len(values), 1))
            probofgettincurrentwindow = sum(np.exp(kde.score_samples(values)))

            essentialregions[chromosome][value['bin']] = {'bin': value['bin'], 'count': sumwindowval, 'essentialpval': probofgettincurrentwindow, 'essentialregions': (1 if probofgettincurrentwindow <= pvalue else 0)}
            if 'regions' not in featurelist[chromosome][value['name_id']]:
                featurelist[chromosome][value['name_id']]['essentialcalculated'] = 0
                featurelist[chromosome][value['name_id']]['essentialpredicted'] = 0
                featurelist[chromosome][value['name_id']]['structure'] = ''
                featurelist[chromosome][value['name_id']]['struc'] = 0
                featurelist[chromosome][value['name_id']]['regions'] = 0
            featurelist[chromosome][value['name_id']]['essentialcalculated'] += (1 if probofgettincurrentwindow <= pvalue else 0)
            featurelist[chromosome][value['name_id']]['regions'] += 1

    return [essentialregions, featurelist]


def discretize5(tnlist):
    discret = {}
    # discret2 = {}
    perfect = {}
    # perfect2 = {}
    for chromosome in sorted(tnlist.keys()):
        dist = []
        discret[chromosome] = []
        # discret2[chromosome] = []
        perfect[chromosome] = []
        # perfect2[chromosome] = []
        for bins in sorted(tnlist[chromosome].keys()):
            dist.append(tnlist[chromosome][bins]['count'])
        thresholds = np.percentile(dist, [1, 25, 50, 75])
        print(thresholds)
        for bins in sorted(tnlist[chromosome].keys()):
            count = tnlist[chromosome][bins]['count']
            if count <= thresholds[0]:
                ret = 0
            elif count <= thresholds[1]:
                ret = 1
            elif count <= thresholds[2]:
                ret = 2
            elif count <= thresholds[3]:
                ret = 3
            else:
                ret = 4
            discret[chromosome].append(ret)
            perfect[chromosome].append(0)
            # discret2[chromosome].append(1 - tnlist[chromosome][bins]['essentialpval'])
            # perfect2[chromosome].append(1)

    return [discret, perfect]  # , discret2, perfect2]


def predict_essential(essentialregions, similation=1000, verbose=False):
    print('>Predict essential regions...')
    # (discta, prefect, discta2, prefect2) = discretize5(essentialregions)
    (discta, prefect) = discretize5(essentialregions)

    for chromosome in sorted(discta.keys()):
        remodel = hmm.GaussianHMM(n_components=5, covariance_type="full", n_iter=similation)
        x = np.array([discta[chromosome]]).T
        y = np.array([prefect[chromosome]]).T
        remodel.fit(x)
        if verbose:
            print(' * Convergence:', str(remodel.monitor_.converged))
        test = remodel.predict(y)
        predicted = remodel.predict(np.array([discta[chromosome]]).T)
        essential = int(sum(test) / len(test))
        for i in range(0, len(discta[chromosome])):
            essentialregions[chromosome][i]['predictedclass'] = predicted[i]
            essentialregions[chromosome][i]['predictedessential'] = bool(predicted[i] == essential)

    return essentialregions


def output_essential(essentialregions, tnwindows, featurelist, prefix, fasta, verbose=False):
    print('>Output results...')
    region = re.compile(r'^(X+O+|O+X+)$')
    regionend = re.compile(r'^O+X+$')
    for chromosome in sorted(essentialregions.keys()):
        for item in essentialregions[chromosome].values():
            if 'predictedessential' in item and item['predictedessential'] is True:
                featurelist[chromosome][tnwindows[chromosome][item['bin']]['name_id']]['essentialpredicted'] += 1
            featurelist[chromosome][tnwindows[chromosome][item['bin']]['name_id']]['structure'] += str('O' if 'predictedessential' in item and item['predictedessential'] is True else 'X')
            featurelist[chromosome][tnwindows[chromosome][item['bin']]['name_id']]['struc'] += int(item['predictedclass'])
    count = {'Neutral': 0, 'Regional': 0, 'Under-represented': 0}

    seqs = {}
    print('>Output summary table...')
    with open(prefix + '.table.tsv', 'w') as handle:
        handle.write('Locus Tag\tLocus\tProtein Name\tTotal TAs\tFraction TAs hit\tTotal reads\tEL-ARTIST\tDisrupted\tGene structure\n')
        for chromosome in sorted(essentialregions.keys()):
            seqs[chromosome] = []
            for entry in sorted(featurelist[chromosome].keys()):
                if not featurelist[chromosome][entry]['name'].startswith('IG_') and 'essentialpredicted' in featurelist[chromosome][entry]:
                    status = 'Neutral'
                    if featurelist[chromosome][entry]['essentialpredicted'] > 0 and featurelist[chromosome][entry]['essentialpredicted'] / featurelist[chromosome][entry]['regions'] > 0.1:
                        if featurelist[chromosome][entry]['essentialpredicted'] == featurelist[chromosome][entry]['regions'] and featurelist[chromosome][entry]['score'] <= 2:
                            status = 'Under-represented'
                        else:
                            m = region.search(featurelist[chromosome][entry]['structure'])
                            r = regionend.search(featurelist[chromosome][entry]['structure'])
                            if m is not None:
                                if r is not None and featurelist[chromosome][entry]['essentialpredicted'] / featurelist[chromosome][entry]['regions'] >= 0.90 and featurelist[chromosome][entry]['score'] <= 4:
                                    status = 'Under-represented'
                                else:
                                    status = 'Regional'
                    handle.write(featurelist[chromosome][entry]['name'] + '\t' + (featurelist[chromosome][entry]['gene'] if 'gene' in featurelist[chromosome][entry] else '-') + '\t' + (featurelist[chromosome][entry]['product'] if 'product' in featurelist[chromosome][entry] else '-') + '\t' + str(featurelist[chromosome][entry]['sites']) + '\t' + (str(round(featurelist[chromosome][entry]['score'] / featurelist[chromosome][entry]['sites'], 2)) if featurelist[chromosome][entry]['sites'] > 0 else '-') + '\t' + str(featurelist[chromosome][entry]['reads']) + '\t' + status + (' (<10 TAs)' if featurelist[chromosome][entry]['sites'] < 10 else '') + '\t' + str(5 * round(20 * (1 - featurelist[chromosome][entry]['essentialpredicted'] / featurelist[chromosome][entry]['regions']), 0)) + '\t' + featurelist[chromosome][entry]['structure'] + '\n')
                    count[status] += 1


    return count


def main(gff, fasta, tn_report, tn_window, prefix, similation=1000, windows=5, pvalue=0.005, feature='gene', site='TA', verbose=False):
    # 1. Count the chromosomes and TN site in every chromosome
    (featurelist, tnlist, chrlist) = parse_genomes(gff, fasta, feature, site, verbose)

    # 2. Count the mapped reads at every TN site in the chromosome
    tnmapped = parse_reports(tn_report, tnlist)

    # 3. Mapping TN insertions into windows
    tnwindows = compact_reports(tn_window, chrlist, tnmapped, featurelist)
    featurelist = update_features(tnmapped, featurelist)

    # 4. Visualizing Sequencing Saturation
    checkseqdepth(tnwindows, prefix, verbose)

    # 5. Normalisation
    # tnwindows = window_average(tnwindows, similation, chrlist)

    # 6. EL-ARTIST:  Essential Loci Analysis ( = smoothing algorithm)
    (essentialregions, featurelist) = slidingwindow(tnwindows, featurelist, similation, windows, pvalue)
    essentialregions = predict_essential(essentialregions, similation, verbose)

    # 7. Output final list
    ret = output_essential(essentialregions, tnwindows, featurelist, prefix, fasta, verbose)
    if verbose:
        print(ret)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', dest='gff', type=str, required=True, help='Path to the gff3/gtf file')
    parser.add_argument('--fasta', dest='fasta', type=str, required=True, help='Path to the fasta file')
    parser.add_argument('--tn', dest='tn_report', type=str, required=True, help='Path to the TN report (gff3/gtf)')
    parser.add_argument('--tn_window', dest='tn_window', type=int, default=10, help='Size of the bin (nucleotids)')
    parser.add_argument('--sl_window', dest='sl_window', type=int, default=5, help='Size of the sliding window the smoothing algorithm')
    parser.add_argument('--simulations', dest='simulations', type=int, default=1000, help='Number of simulation to run')
    parser.add_argument('--pvalue', dest='pvalue', type=float, default=0.005, help='P-value threshold for "essential calculated"')
    parser.add_argument('--feature', dest='feature', type=str, default='CDS', help='Feature to retain: e.g. gene or CDS')

    parser.add_argument('--prefix', dest='prefix', type=str, default='output', help='Output prefix')
    parser.add_argument('--verbose', dest='verbose', action='store_true', help='Become very chatty')

    args = parser.parse_args()

    if args.gff is not None and os.path.exists(args.gff) and args.fasta is not None and os.path.exists(args.fasta) and args.tn_report is not None and os.path.exists(args.tn_report):
        main(args.gff, args.fasta, args.tn_report, int(args.tn_window), args.prefix, args.simulations, args.sl_window, args.pvalue, args.feature, 'TA', args.verbose)
    else:
        parser.print_help()
        exit(1)
