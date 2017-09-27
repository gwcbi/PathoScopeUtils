#! /usr/bin/env python
import sys
import argparse

samflags = {
    'PAIRED': 0x1,
    'PROPER_PAIR': 0x2,
    'UNMAP': 0x4,
    'MUNMAP': 0x8,
    'REVERSE': 0x10,
    'MREVERSE': 0x20,
    'READ1': 0x40,
    'READ2': 0x80,
    'SECONDARY': 0x100,
    'QCFAIL': 0x200,
    'DUP': 0x400,
    'SUPPLEMENTARY': 0x800
}

def sam_generator(fh):
    samlines = (l.strip('\n').split('\t') for l in fh)
    header = []
    for sl in samlines:
        if sl[0].startswith('@'):
            header.append(sl)
        else:
            yield header
            break
    
    alns = [sl,]
    for sl in samlines:
        if sl[0] == alns[0][0]:
            alns.append(sl)
        else:
            yield alns
            alns = [sl, ]
    
    # Final set:
    if alns:
        yield alns

# Functions for converting SAM tags
conv_funcs = {
    'A': lambda x: str(x)[0],
    'i': int,
    'f': float,
    'Z': str,
    'H': lambda x: [int(x[i:i+2],16) for i in range(0,len(x),2)],
    'B': lambda x: map(float if x[0]=='f' else int, x[1:].split(',')),
}

def parse_sam_tags(samline):
    """ Return dictionary with all tags """
    return {t[:2] : conv_funcs[t[3]](t[5:]) for t in samline[12:]}

def get_sam_tag(samline, tag):
    """ Return value for tag """
    for t in samline[12:]:
        if t[:2] == tag:
            return conv_funcs[t[3]](t[5:])
    return None


def main(args):
    sg = sam_generator(args.infile)
    header = sg.next()
    for hl in header:
        print >>args.outfile, '\t'.join(hl)
    sec_to_pri = 0
    for alns in sg:
        if len(alns) == 1:
            aln = alns[0]
            if int(alns[0][1]) & samflags['SECONDARY']:
                alns[0][1] = '%d' % (int(alns[0][1]) ^ samflags['SECONDARY'])
                sec_to_pri +=  1
            print >>args.outfile, '\t'.join(alns[0])
        else:
            # First sort by reference name
            salns = sorted(alns, key=lambda x: x[2])
            # Then sort by alignment score, largest to smallest
            salns.sort(key=lambda x: get_sam_tag(x, 'AS'), reverse=True)
            # Change first alignment to primary
            if int(salns[0][1]) & samflags['SECONDARY']:
                salns[0][1] = '%d' % (int(salns[0][1]) ^ samflags['SECONDARY'])
            
            # Change other alignments to secondary
            for i in range(1,len(salns)):
                salns[i][1] = '%d' % (int(salns[i][1]) | samflags['SECONDARY'])
            
            print >>args.outfile, '\n'.join('\t'.join(a) for a in salns)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix SAM flags for unpaired reads')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), 
                        default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)
    main(parser.parse_args())
