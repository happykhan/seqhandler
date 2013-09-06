SeqHandler
========

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

CHANGE LOG
----------
2013-04-16 Nabil-Fareed Alikhan <n.alikhan@uq.edu.au>
    * Version 0.3 
    * Initial build
2013-04-17 Nabil-Fareed Alikhan <n.alikhan@uq.edu.au> 
    * Version 0.5 
    * Reworked GBKSplit to SeqHandler
    * Created sub modules: split, merge & convert files
    * Explicit control of Input/Output for modules
2013-09-05 Nabil-Fareed Alikhan <n.alikhan@uq.edu.au>-
    * Changed fasta header handling for Prokka input in merge
    * Added header override flags for merge
2013-09-05 Mitchell Stanon-Cook <m.stantoncook@uq.edu.au>-
    * Made into an installable package
    * Installs a script (SeqHandler) system wide
    * Small improvements in terms of using __init__ as a meta container
2013-09-06 Mitchell Stanon-Cook <m.stantoncook@uq.edu.au>-
    * Added option to convert to/from gff


LICENCE
=======
Nabil-Fareed Alikhan <n.alikhan@uq.edu.au>. (C) 2012-2013.

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

