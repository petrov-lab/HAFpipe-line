#!/usr/bin/env python

import sys, os

# Get the input snp table.
if len(sys.argv) != 2:
    sys.exit("Need snptable file name")
snptable_file = sys.argv[1]

# Get the two derived file names, allele counts for input, and numeric for output.
# Allele count file ac_file must be in the same row-order as snptable_file,
# and have a column for identity of the alt allele.
ac_file = snptable_file + ".alleleCts"
numeric = snptable_file + ".numeric"

# Open all files.
snptable_in = open(snptable_file, 'r')
ac_in       = open(ac_file, 'r')
numeric_out = open(numeric, "w")

# Get the header lines, and check them for what we expect.
snptable_header = snptable_in.readline().rstrip().split(',')
ac_header_line  = ac_in.readline()
if snptable_header[1] != "Ref":
    raise Exception("snptable file " + snptable_file + " does not start with the expected header")
if ac_header_line.rstrip() != "pos,ref,alt,nCt,refCt,altCt,hetCt":
    raise Exception("alleleCts file " + ac_file + " does not start with the expected header")

# Write the header of our output file.
# We just skip the "Ref", the rest is the same as in the original table.
chrom = snptable_header[0]
numeric_out.write( chrom + ',' + ','.join( snptable_header[2:] ) + '\n' )

# Now iterate both files in parallel.
for snp_line, ac_line in zip( snptable_in, ac_in ):
    snps = snp_line.rstrip().split(',')
    acs  = ac_line.rstrip().split(',')

    # We except:
    #   - snps[0] and acs[0]: position
    #   - snps[1] and acs[1]: reference base
    #   - acs[2]: alternative base
    #   - snps[2:]: bases as found in the samples

    # Check that we are in the same position in the genome in both tables,
    # and that the ref bases are identical.
    if snps[0] != acs[0]:
        raise Exception(
            "snptable file " + snptable_file + " and its .alleleCts file are out of sync; " +
            "the former is at " + str( snps[0] ) + " while the latter is at " + str( acs[0] )
        )
    if snps[1] != acs[1]:
        raise Exception(
            "snptable file " + snptable_file + " and its .alleleCts file are out of sync; " +
            "the former has ref " + str( snps[1] ) + " while the latter has ref " + str( acs[1] ) +
            " at position " + str( snps[0] )
        )

    # Write position to the output file.
    numeric_out.write( snps[0] )

    # Compute and write the numeric values for each snp of each sample.
    # Code: ref=0; alt=1; het=.5; missing=-1
    for snp in snps[2:]:
        # Formula to convert this elegantly to these values is taken from the original R script.
        # We are not using this, as the below explicit version is 4 times faster...
        # is_r = int( snp == acs[1] ) # compare to ref allele
        # is_a = int( snp == acs[2] ) # compare to alt allele
        # is_n = int( snp == 'N' )
        # val = 0.5 - 0.5 * is_r + 0.5 * is_a - 1.5 * is_n
        # numeric_out.write( ',' + '%g'%( val ))

        if snp == acs[1]:
            numeric_out.write( ',0' )
        elif snp == acs[2]:
            numeric_out.write( ',1' )
        elif snp == 'N':
            numeric_out.write( ',-1' )
        else:
            numeric_out.write( ',0.5' )

    # Finish the line
    numeric_out.write( '\n' )

# Finished
snptable_in.close()
ac_in.close()
numeric_out.close()
