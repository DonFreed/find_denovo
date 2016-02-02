find_denovo
===========

A simple tool for identifying *de novo* variants from a VCF/BCF file.

### Installation
```
git clone git://github.com/DonFreed/find_denovo.git
make HTSDIR=/path/to/htslib/
```

### Usage
mark_denovo [options] [-O <v|z|b|u>] [-i <in.vcf|in.vcf.gz|in.bcf>] -p in.ped [-o <out.vcf|out.vcf.gz|out.bcf>]

Options:
           -c INT          Minimum number of reads in all trio members [20]
           -a INT          Minimum number of reads supporting the alternate allele [3]
           -t INT          Minimum phred-scaled confidence for the parental genotype [20]
           -s INT          Minimum phred-scaled confidence for the child's genotype [20]
           -l INT          zlib compression level for the output VCF/BCF [7]
           -i FILE         The input file [stdin]
           -p FILE         The input PED file
           -o FILE         The output file [stdout]
           -d FILE         Output abbreviated information on the identified de novo variants (recommended for large inputs)
           -O <v|z|b|u>    v: VCF, z: bgzip compressed VCF, b: BCF, u: uncompressed BCF [v]
           -n INT          The maximum number of non-sibling individuals with the allele for a de novo call [2]
           -x INT          number of extra compression/decompression threads [0]
           -h              Print this help information
