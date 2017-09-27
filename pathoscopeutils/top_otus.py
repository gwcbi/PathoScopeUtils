#! /usr/bin/env python

import sys
import argparse
import re

def get_top_otus(f, cutoff=0.01, with_tax=False):
    ''' Get the top OTUs in a PathoScope report '''
    rows = [l.strip('\n').split('\t') for l in open(f, 'rU')]
    h1,h2 = rows[:2]
    rows = rows[2:]
    tot = sum(float(r[3]) for r in rows)
    top_rows = [r for r in rows if float(r[3]) > (cutoff * tot) ]
    # Strip off ti
    tis = {}
    for r in top_rows:
        m = re.search('ti\|(\d+)\|?', r[0])
        ti = m.group(1) if m else r[0]
        if not with_tax:
            tis[ti] = True
        else:
            taxname = ''
            idx = -1
            while not taxname:
                taxname = r[idx]
                idx += -1
            tis[ti] = taxname
    
    return tis

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Return top OTUs for PathoScope report(s)')
    parser.add_argument('--with_tax', action='store_true')
    parser.add_argument('--pct', type=float, default=0.01)
    parser.add_argument('reports', nargs='+')
    args = parser.parse_args()
    
    top_otus = {}
    for r in args.reports:
        top_otus.update(get_top_otus(r, args.pct, args.with_tax))
    
    if not args.with_tax:
        print >>sys.stdout, '\n'.join(sorted(top_otus.keys()))
    else:
        for k in sorted(top_otus.keys()):
            print >>sys.stdout, '%s\t%s' % (k, top_otus[k])
