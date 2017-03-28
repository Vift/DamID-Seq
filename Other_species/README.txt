If you wish to use DamID-seq pipeline with any other (not Drosophila melanogaster) species, you should use or make this species-specific "All_GATC_list.txt" file.
Precomputed files for human, mouse, C.elegans and D.rerio can be found in the corresponding folders. They must be uncompressed and placed in the same 
folder with the scripts or executables (with replacement of an old Drosophila-specific "All_GATC_list.txt" file).
You can make "All_GATC_list.txt" for any other species with the Perl script "GATC_coordinates_for_genome.pl". It requires single-file genome in fasta format in the 
same folder with the script. File must have extension ".fa" or ".fa.gz". In the case of compressed genome file make sure that the archive doesn't contain any folders,
only one .fa file. Alternatively you can write an inquiry for "All_GATC_list.txt" for a custom genome to "Vift@mcb.nsc.ru", I'll try to make it and add it to this depository.