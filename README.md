0;136;0c# ConsDB
A Python tool for interfacing with large variant databases and performing consensus genome operations.

Written by Benjamin Kaminow

# Usage
The main ConsDB script can be run as follows:

    cd consdb
    python ConsDB.py <run mode> [arguments]
    
Where `<run mode>` is one of the following:
*  `Parse` - Parse database files into ConsDB files
*  `Filter` - Filter a VCF file to remove major alleles
*  `Merge` - Merge multiple ConsDB files
*  `Cons` - Create a consensus VCF file
*  `FA` - Create a consensus FASTA file

Arguments for each run mode can be viewed by passing `-h` as the sole argument.

# Back-End
ConsDB utilizes an object oriented back-end, which can easily be incorporated into pipeline/scripts. 
Documentation for this back-end is found in the `docs` subfolder.
