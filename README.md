# *T-Aligner*
------------

U-insertion/deletion editing study tool.

Current release

*T-Aligner v.3.1.mt*

Changes:

> - multithreaded on most steps
> - fixed bug with seed mismatch detection
> - improved stability

Binary files updated!

## Running *T-Aligner*
------------

 *T-Aligner* is command-line tool to analyze RNA editing in
Trypanosomatids. To run T-Aligner you will need reference
file in `.fasta` format and reads file in `.fastq` format.
Reference should contain single cryptogene sequence written
only using A/C/G/T letters alphabet.

Example:

```
>CYb
TAACTATTTAAAATAATAATATATGTGATATATGAAT
TATTAAAATAAAAAGCGGAGAAAAGAGGAAGGGTCTT
CTAATGTCAG
```

Reference sequence should better be flanked with 50-100 bp
regions.

Typical launch command to map reads on reference sequence
and get ORFs `fasta` files, alignments and pictures is:

>./T-Aligner-3 --in_ref path/to/reference.fasta --in_lib path/to/reads.fastq --ss 10 --sl 10 --j 50000 --t 2 --out output_prefix

This is a minimal command line with mandatory parameters. `--ss` and `--sl` specify seed parameters and it is highly
recommended to set ss=sl. Usually you DO NOT need to change seed step and seed length!
`--j` and `--t` are technical parameters that specify number of reads per job and number of threads.
You can leave `--j` unchanged, but specify the number of threads suitable for you machine.
Output prefix (`--out`) will be used for naming output files.

There are also plenty of other optional parameters you can set:

`--dump_mapped_reads` - save fraction of mapped reads in separate .fastq file. Usefull when you will rerun *T-Aligner* on the same dataset.

`--draw_orf [int ID]` - draw ORF with ID (default: 0 - main ORF, longest)

`--mf [float]` - minimal mapped fraction of read's length (default 0.8)

`--mr [int]` - minimal read length (default 16, should be at least of seed length)

`--start_with_atg` - allow valid ORF to start only from ATG (default ATG/ATA)

`--orf_length_min [int ID]` - set minimal ORF length

`--rg_overlap_min` - minimal overlap length when ORF assembly graph construction

`--draw_orf_rgb [int ID],[int R],[int G],[int B]` - draw ORF with ID and RGB color set with 3 integers 0-255 [int R],[int G],[int B]

`--draw_orfs_btlc [int X],[int Y],[int K]` - draw K longest ORFs passing through an editing matrix point with T-less coordinates X,Y.


## How to build
------------

You will need Qt library to build *T-Aligner*.
Download and install Qt, make sure that qmake
is in your environment's PATH.
Download T-Aligner's code and do:

>qmake -project

>qmake

>make

Some versions of GCC/C++ complier need to
specify `-std=c++11 -stdlib=libc++` flags in Makefile.
To do this please open Makefile in any text editor and
add '-std=c++11 -stdlib=libc++' to CFLAGS and CXXFLAGS.


>CXXFLAGS=-std=c++11 -stdlib=libc++ ...

>CFLAGS=-std=c++11 -stdlib=libc++ ...


Usually no other corrections needed. If you still
experience any problems with building *T-Aligner*
please use pre-built executable files available
in /bin directory specific to your OS.

