#!/usr/bin/env python3
"""
Extract strains carrying each mutation in a VCF file
"""
import argparse
import fileinput
import vcf

CHR_MAP = {'chromosome1': 'chrI', 'chromosome2': 'chrII', 'chromosome3': 'chrIII',
           'chromosome4': 'chrIV', 'chromosome5': 'chrV', 'chromosome6': 'chrVI',
           'chromosome7': 'chrVII', 'chromosome8': 'chrVIII', 'chromosome9': 'chrIX',
           'chromosome10': 'chrX', 'chromosome11': 'chrXI', 'chromosome12': 'chrXII',
           'chromosome13': 'chrXIII', 'chromosome14': 'chrXIV', 'chromosome15': 'chrXV',
           'chromosome16': 'chrXVI'}

def main(args):
    """Main script"""
    # Import strain info
    if args.meta:
        meta = {}
        with fileinput.input(args.meta) as meta_file:
            for line in meta_file:
                if meta_file.isfirstline():
                    header = line.strip().split('\t')
                    standard = header.index('Standardized name')
                    isolate = header.index('Isolate name')

                else:
                    line = line.strip().split('\t')
                    meta[line[standard]] = line[isolate]

    vcf_file = vcf.Reader(filename=args.vcf)

    # Print headers
    if args.meta:
        print('chromosome', 'position', 'ref', 'alt', *[meta[i] for i in vcf_file.samples], sep='\t')

    else:
        print('chromosome', 'position', 'ref', 'alt', *vcf_file.samples, sep='\t')

    # Iterate over loci
    for site in vcf_file:
        freqs = site.INFO['AF']
        for i in range(len(site.ALT)):
            if not args.filter or float(freqs[i]) < args.filter:
                gens = [call.data.GT.count(str(i + 1)) for call in site.samples]
                print(CHR_MAP[site.CHROM], site.POS, site.REF, site.ALT[i], *gens, sep='\t')

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('vcf', metavar='V', help="Input VCF file")

    parser.add_argument('--meta', '-m', default='',
                        help="File giving strain details to convert strain names to \
                              standardized names used in VCF")

    parser.add_argument('--filter', '-f', default=0.0, type=float,
                        help="Filter to genotypes at frequencies lower than the given level")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
