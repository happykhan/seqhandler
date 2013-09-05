#!/usr/bin/env python
"""
# Created: Tue, 16 Apr 2013 14:49:03 +1000

SeqHandler provides low level operations on sequence files (FASTA, GenBank or 
EMBL files).

Requires Biopython, see: http://biopython.org/wiki/Main_Page 
& Python 2.7.3 or similar

SPLIT MODULE
============
Use this Module if you want to split a gbk/embl file on a particular feature.
e.g. source, fasta_record

An example would be :
python SeqHandler.py split test.embl splitTest -f source -i embl

Where:
    * test.embl is the input file you want to split (gbk or embl). 
    * -i flag, specifies the input file type (genbank/embl)
      If none is specified the script will assume input is a  Genbank file
    * splitTest is the output directory for all the resulting files.
    * -o flag, specifies the output file type (genbank/embl/fasta)
    * -f flag, denotes the feature you want to split on (In this case 
      'source'). If none is specified the script will split on 'fasta_record'

N.B The script will try to rename the split files using the 'note' attached 
to the feature you're trying to split on. 

If no output format is specified (-o) flag, the file format for split files 
will be the same as the input. i.e embl input, get embl output.

MERGE MODULE
============
Use this module if you want to merge a genbank, embl or fasta records. 

An example would be: python SeqHandler.py merge test.gbk -i genbank test2.fna -o fasta

Where:
    * test.gbk is the input file you want to merge (gbk/embl/fasta). 
    * test2.fna is the output file.
    * -i flag, specifies the input file type (genbank/embl/fasta)
      If none is specified the script will assume input is a  Genbank file
    * -o flag, specifies the output file type (genbank/embl/fasta)

N.B Input only takes one file. If you want to merge multiple seperate files
together you'll need to first concatenate them together. 

If no output format is specified (-o) flag, the file format for split files 
will be the same as the input. i.e embl input, get embl output.


CONVERT MODULE
==============
Use this module if you want to convert a sequence file. 

An example would be: 
python SeqHandler.py convert test.gbk test2.fna -o fasta -i genbank

Where:
    * test.gbk is the input file you want to merge
    * test2.fna is the output file.
    * -i flag, specifies the input file type 
      If none is specified the script will assume input is a Genbank file
    * -o flag, specifies the output file type
      If none is specified the script will assume input is a fasta file

See USAGE (python SeqHandler.py convert -h) for full list of supported files.

### CHANGE LOG ### 
2013-04-16 Nabil-Fareed Alikhan <n.alikhan@uq.edu.au>
    * Version 0.3 
    * Initial build
2013-04-17 Nabil-Fareed Alikhan <n.alikhan@uq.edu.au> 
    * Version 0.5 
    * Reworked GBKSplit to SeqHandler
    * Created sub modules: split, merge & convert files
    * Explicit control of Input/Output for modules
2013-09-05 Nabil-Fareed Alikhan <n.alikhan@uq.edu.au> 
    * Changed fasta header handling for Prokka input in merge
    * Added header override flags for merge
2013-09-05 Mitchell Stanon-Cook <m.stantoncook@uq.edu.au>-
    * Made into an installable package
    * Installs a script (SeqHandler) system wide
    * Small improvements in terms of using __init__ as a meta container

"""
import SeqHandler.__init__ as meta
import sys, os, traceback, argparse
import time
from Bio import SeqIO

epi = "Licence: %s by %s <%s>" % (meta.__license__, 
                                  meta.__author__,
                                  meta.__author_email__)
__doc__ = " %s v%s - %s (%s)" % (meta.__title__, 
                                 meta.__version__, 
                                 meta.__description__, 
                                 meta.__url__)

def convertMod(args):
    try:
        count = SeqIO.convert(args.input, args.inFormat, args.output, args.outFormat )
        if count == 0:
            sys.err.write('ERROR: No records converted. Possibly wrong input filetype\n')
        else: 
            if args.verbose: sys.err.write("Converted %i records\n" %count)
    except ValueError:
        sys.err.write("ERROR: ValueError, Could not convert file\n")

def mergeMod(args):

    filetype = args.inFormat
    # Load file as SeqRecord
    int_handle = open(args.input, "r")
    recs = list(SeqIO.parse(int_handle, filetype))
    # For each SeqRecord, I.e. complete gbk annotation obj in file
    fgbk = recs[0]
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    d = SeqFeature(FeatureLocation(0, len(fgbk) ), type="fasta_record",\
                strand=1)
    d.qualifiers["note"] = recs[0].name 
    fgbk.features.append(d)
    for l in recs[1:]:
        d = SeqFeature(FeatureLocation(len(fgbk), len(fgbk) + len(l)), type="fasta_record",\
                strand=1)
        d.qualifiers["note"] = l.name 
        fgbk.features.append(d)
        fgbk += l
    fgbk.name = recs[0].name
    fgbk.description = recs[0].description
    fgbk.annotations = recs[0].annotations
    if args.accession != None: 
        fgbk.name = args.accession
    if args.ver != None:
        fgbk.id = fgbk.name +'.' + args.ver
    for f in fgbk.features:
        if f.type == 'source':
            fgbk.features.remove(f)
    d = SeqFeature(FeatureLocation(0, len(fgbk)), type="source", strand=1)
    fgbk.features.insert(0,d)
    outtype = filetype 
    if args.outFormat != None:
        outtype = args.outFormat
    out_handle = open( args.output,"w")
    SeqIO.write(fgbk, out_handle, outtype)

def splitMod(args):

    # Input Filetype, check if gbk or embl
    filetype = args.inFormat
    # Load file as SeqRecord
    int_handle = open(args.input, "r")
    Recs = SeqIO.parse(int_handle, filetype)
    if not os.path.exists(args.outputDir):
        os.mkdir(args.outputDir)
    # For each SeqRecord, I.e. complete gbk annotation obj in file
    count = 0 
    for rec in Recs:
        for feat in rec.features:
            # Split file everytime we see user-defined feature 
            # (fasta_record if not specified, could be source if user 
            # Uses ' -f source ' 
            if feat.type == args.feature:
                count += 1 
                subname = '' 
                # Append note/annotation for feature to filename 
                #TODO: Clean up note field is pretty basic
                if feat.qualifiers.has_key('note'):
                    subname = feat.qualifiers['note'][0].replace(' ','-')\
                            .replace(';','')
                # TODO: ID/Accession  would probably need to have contig 
                # number suffixed
                Finalname = subname +  str(os.path.basename(args.input))
                outtype = filetype
                if args.outFormat != None:
                    outtype = args.outFormat
                outSuffix = '.' + outtype 
                if outtype == 'genbank': outSuffix = '.gbk'
                if not Finalname.endswith(outSuffix): Finalname += outSuffix
                out_handle = open( os.path.join(args.outputDir, Finalname ),\
                        "w")
                # Create new SeqRecord between location of split feature
                finalGbk = rec[feat.location.start:feat.location.end]
                # Copy annotations-Header from old file to new file
                finalGbk.annotations = rec.annotations
                # Write Output file to specified dir 
                SeqIO.write(finalGbk, out_handle, outtype)
    if count == 0:
        sys.write.err('No file generated; Wrong feature specified? ' + 
        'see -f flag\n')

if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = argparse.ArgumentParser(description=__doc__ ,epilog=epi)
        parser.add_argument ('-v', '--verbose', action='store_true', \
                default=False, help='verbose output')
        parser.add_argument('--version', action='version', version='%(prog)s '\
                + meta.__version__)
        subparsers = parser.add_subparsers(help='modules')
        split_parser = subparsers.add_parser('split', help='Splits sequence files')
        split_parser.add_argument('-f','--feature',action='store',default=\
                'fasta_record',\
                help='Annotation feature to split on [def: "fasta_record"]')
        split_parser.add_argument ('input', action='store', \
                help='Input annotation file to split')
        split_parser.add_argument('-i','--inFormat',action='store',\
                choices=('embl', 'genbank'),\
                default='genbank', help='Format of input file [def: genbank]')
        split_parser.add_argument('-o','--outFormat',action='store',\
                choices=('embl', 'fasta', 'genbank'),\
                help='Format of output file')
        split_parser.add_argument ('outputDir', action='store', \
                help='Output Directory')
        merge_parser = subparsers.add_parser('merge', help='Merges sequence files')
        merge_parser.add_argument ('input', action='store', \
                help='Input annotation file')
        merge_parser.add_argument ('output', action='store', \
                help='Output File')
        merge_parser.add_argument('-i','--inFormat',action='store',\
                choices=('embl', 'fasta', 'genbank'),\
                default='genbank', help='Format of input file [def: genbank]')
        merge_parser.add_argument('-a','--accession',action='store',\
                help='User defined accession no/locus')
        merge_parser.add_argument('-e','--ver',action='store',\
                help='User defined version no')
        merge_parser.add_argument('-o','--outFormat',action='store',\
                choices=('embl', 'fasta', 'genbank'),\
                help='Format of output file')
        convert_parser = subparsers.add_parser('convert', help='Converts sequence files')
        convert_parser.add_argument ('input', action='store', \
                help='Input annotation file')
        convert_parser.add_argument ('output', action='store', \
                help='Output File')
        convert_parser.add_argument('-i','--inFormat',action='store',\
                choices=('abi', 'ace', 'clustal', 'embl', 'fasta', 'fastq',\
                'fastq-solexa','fastq-illumina','genbank','ig','imgt','nexus',\
                'phd','phylip','pir','seqxml','sff','stockholm','swiss','tab'\
                'qual','uniprot-xml'),
                default='genbank', help='Format of input file [def: genbank]')
        convert_parser.add_argument('-o','--outFormat',action='store',\
                choices=('clustal', 'embl', 'fasta', 'fastq','fastq-solexa',\
                'fastq-illumina','genbank','imgt','nexus', 'phd','phylip',\
                'seqxml','sff','stockholm','tab', 'qual'),\
                default='fasta', help='Format of output file [def: fasta]')
        split_parser.set_defaults(func=splitMod)
        convert_parser.set_defaults(func=convertMod)
        merge_parser.set_defaults(func=mergeMod)
        args = parser.parse_args()
        args.func(args)
        if args.verbose: print "Executing @ " + time.asctime()
        if args.verbose: print "Ended @ " + time.asctime()
        if args.verbose: print 'total time in minutes:',
        if args.verbose: print (time.time() - start_time) / 60.0
        sys.exit(0)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)

