#! /usr/bin/env python

import sys
import os
import argparse
import re
from subprocess import Popen, PIPE, check_call, check_output, CalledProcessError
import subprocess
from glob import glob
from collections import defaultdict

def simplify_filenames(filenames):
    tmp = unique_substring(filenames, '/')
    tmp = unique_substring(tmp, '.')
    tmp = unique_substring(tmp, '_')
    return tmp

def unique_substring(slist, sep):
    """ Return a substring that is unique for all strings """
    z = [s.split(sep) for s in slist]
    for i in range(len(z[0])-1, -1, -1):
        f = [_[i] for _ in z]
        if len(slist) == len(set(f)):
            return f
    return filenames

def main(args):
    bamfiles = sorted(args.bamfile)
    bamnames = simplify_filenames(bamfiles)
    
    print >>sys.stderr, '\n***Sorted BAM files:'
    print >>sys.stderr, '\t%s%s' % ('Name'.ljust(30), 'BAM file')
    for b,n in zip(bamfiles, bamnames):
        print >>sys.stderr, '\t%s%s' % (n.ljust(30), b)

    if args.taxfile:
        print >>sys.stderr, '\n*** Parsing taxonomy: %s' % args.taxfile.name
        top_otus = [l.strip('\n').split('\t') for l in args.taxfile]
        if all(len(l)==2 for l in top_otus):
            top_otus = dict(top_otus)
        else:
            top_otus = {l[0]:l[0] for l in top_otus}
        
        print >>sys.stderr, '\n'.join('%s\t%s' % (k, top_otus[k]) for k in sorted(top_otus.keys()))
    else:
        if args.taxname:
            top_otus = {args.taxid: args.taxname}
        else:
            top_otus = {args.taxid: args.taxid}

    #--- Get original header
    print >>sys.stderr, '\n*** Parsing header...'    

    p1 = Popen('samtools view -H %s' % bamfiles[0], stderr=PIPE, stdout=PIPE, shell=True)
    o, e = p1.communicate()
    header = [l.split('\t') for l in o.strip('\n').split('\n')]
    h1 = header[0]
    
    by_ti = defaultdict(list)
    for hl in header:
        if hl[0] == '@SQ':
            sq_dict = dict(f.split(':') for f in hl[1:])
            m = re.search('ti\|(\d+)\|?', sq_dict['SN'])
            assert m is not None, 'Could not find ti: %s' % '\t'.join(hl)
            if m.group(1) in top_otus:
                by_ti[m.group(1)].append(hl)
    
    #--- Create regions and header files for each OTU
    print >>sys.stderr, '\n*** Creating OTU files...'
    
    for ti, hlines in by_ti.iteritems():
        longname = re.sub(r'\W+', '_', top_otus[ti])
        dirname = os.path.join(args.outdir, longname)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        regfile = os.path.join(dirname, 'regions.bed')
        headfile = os.path.join(dirname, 'header.sam')    
        with open(regfile, 'w') as outR, open(headfile, 'w') as outH:
            print >>outH, '\t'.join(h1)
            for hl in hlines:
                print >>outH, '\t'.join(hl)        
                d = dict(f.split(':') for f in hl[1:])
                print >>outR, '%s\t0\t%s' % (d['SN'], d['LN'])
        
        xbams = []
        for bf,samp in zip(bamfiles, bamnames):
            print >>sys.stderr, '\n***Extracting %s from %s...' % (longname, bf)
            t1 = '%s/tmp1.%s.sam' % (dirname, samp)
            t2 = '%s/tmp2.%s.bam' % (dirname, samp)        
            final_bam = '%s/%s.bam' % (dirname, samp)
            
            # Subset the BAM file
            rgcmd = "-r 'ID:%s' -r 'LB:%s' -r 'PL:SeqLL2' -r 'SM:%s' -r 'PU:%s'" % (samp,samp,samp,samp)
            cmds = [
                'cat %s > %s' % (headfile, t1),
                'samtools view -L %s %s >> %s' % (regfile, bf, t1),
                'samtools view -ub %s | samtools sort > %s' % (t1, t2),
                'picard AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=SeqLL2 RGSM=%s RGPU=%s' % (t2, final_bam, samp, samp, samp, samp),
                'samtools index %s' % final_bam,
            ]
            for cmd in cmds:
                print >>sys.stderr, cmd
                z = check_call(cmd, shell=True)
            
            os.remove(t1)
            os.remove(t2)
            xbams.append(final_bam)
        
        print >>sys.stderr, '\n***Check output in %s' % dirname
        print >>sys.stderr, '\n'.join('\t%s' % xbf for xbf in xbams)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Split BAM files by OTU',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--outdir', default='otu_analysis',
                        help="Output directory")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--taxid', 
                        help="Taxonomy ID to be extracted")
    parser.add_argument('--taxname',
                        help="Name for taxon, used with --taxid")                        
    group.add_argument('--taxfile', type=argparse.FileType('r'),
                        help='''File containing taxonomy IDs to be extracted. If
                                a tab-delimited file, the first column is expected to be
                                the taxonomy IDs and the second column contains a name
                                for the OTU.''')
    parser.add_argument('bamfile', nargs='+',
                        help="Sorted bamfile(s) for extraction")
    args = parser.parse_args()
    
    # Check for executables
    try:
        o,e = Popen('samtools', stderr=PIPE).communicate()
        m = re.search('Version:\s+([\d\.]+)', e)
        if int(m.group(1).split('.')[0]) < 1:
            print >>sys.stderr, "ERROR: samtools must be version 1 or greater"
            print >>sys.stderr, e
            sys.exit()
    except OSError:
        print >>sys.stderr, "ERROR: samtools executable could not be found"    
        sys.exit()
    
    try:
        o,e = Popen('picard', stderr=PIPE).communicate()
    except OSError:
        print >>sys.stderr, "ERROR: picard executable could not be found"    
        sys.exit()
    
    # Run it!
    main(args)
