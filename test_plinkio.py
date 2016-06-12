from plinkio import plinkfile

plink_file = plinkfile.open( "simulated.test" )
if not plink_file.one_locus_per_row( ):
     print( "This script requires that snps are rows and samples columns." )
     exit( 1 )

sample_list = plink_file.get_samples( )
locus_list = plink_file.get_loci( )

for locus, row in zip( locus_list, plink_file ):
    #for sample, genotype in zip( sample_list, row ):
    #    print( "Individual {0} has genotype {1} for snp {2}.".format( sample.iid, genotype, locus.name ) )
    print locus.name,row
    #print "\n"