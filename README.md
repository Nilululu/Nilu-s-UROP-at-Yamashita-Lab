Created on Thu Jun 25 16:45:26 2026


Hi, in this README I will give a overall review of my code structure, all the modules, 
and the order they work. I will also indicate if I need help or feedback on anything specific


ncbi_api.py:
    
    Summery: an script to download all refseq annotated genome folders from NCBI using a tsv file with 
    the accesion ids of genome folders and a API link

    files and folders created by running it: 
        ncbi_data_directory_0.1, genomic_directory.csv, error_dowload_log.txt, id_set.txt
    
    unzips and organize all successfully downloaded folders in ncbi_data_directory_0.1/
    
    logs warnings and errors in error_download_log.txt
    
    finds the location of the gtf file in the downlaoded folder and saves it in genomic_directory.csv
    
    keeps track of accession ids of all successfully downloaded genomic folders and save them in id_set.txt

ncbi_api_2.py:
   is the same as ncbi_api.py except that it downloads independent annotations!
   files and folders created by running it: 
       ncbi_data_directory_0.2, genomic_directory_gca.csv, error_dowload_log2.txt, id_set_gca.txt
    
extract.py:
    
    Summery: a module to take a gtf file as input and extract the information of interest about all genes in a genome 
    uses a dictionary data structure
    
    Data Structure:
        {gene_id: {"position": (start, end)", "transcrips":{}, "introns":set(), "exons":set(), "strand": , "products":set(), "chr":}}
        
            the "transcript" key has the follwing structure to its associated value
            {transcript_id: {"exons":[(start, end), ...] , "introns":[(start, end), ...]}, "product": }
    
    It works in two steps!
    To use it you should call extract_id_and_genes(gtf_file), this will give a dictionary with the data structure mentioned above
    Then call compute_introns on the said dictioanry, this will use the information about exons to compute and store information about introns
        
        
taxonomy.py:
    a module to take a taaxonomy number and use taxonomy dataset to return a taxonomy dictionary
    It needs names.dmp and nodes.dmp to be in the working directory / have specified locations
    
    Functions:
        generate_taxonomy_dict:
            uses nodes.dmp to make the following dictioanry
            taxonomy_dict[tax_id] = (parent_tax_id, rank)
            
        generate_tax_to_names:
            uses the names.dmp file to make the folowing dictionary 
            tax_to_name[tax_id]= name
        
        find_taxonomy:
            uses a specific tax_num and the dictionaries created by the previous two functions 
            to create a toxonomy lineage dictionary 
            
metadata.py:
    find_g_intron_info: (! nort currently used since we are still trying to define what makes an intron a giant intorn)
        gets a dictionary created by extract module and uses a threshold to identify giant introns 
        and their realtivve start and end position compared to the gene they are on
        
     get_genome_metadata:
         uses the gtf file location to find the assembly report jsonl file 
         makes a josn dictionary out of the assemblt report 
         seraches for keys of interest (tax_id, total_sequence_length, assembly_level, assembly_type, 
         numChr, num_scaffolds, num_contigs, scaffold_n50, contig_n50, gc_percent) 
    find_key_multiple:
        a helper functiont for get_genome_metadata to search for keys of interest


plot_creator: OUTDATED (! right now not used, useful if you want to make the same plot for many genomes 
               without creating a plot object in your code again and again)
    Functions:
        create_histogram
        create_scatter (3d, incorporates colors as a 3rd axis)
        create_2d_scatter
         
    
test_part_1.py:
    
    Summery: uses genmic_directory.csv and extract, taxonomy, and metadata modules tp extract the 
    following information and write it in a text file :
        genome_id, name, kingdom, tax_id, total_sequence_length, assembly_level, assembly_type, 
        numChr, num_scaffolds, num_contigs, scaffold_n50, contig_n50, gc_percent, max_intron, 
        min_intron, mean_intron,median_intron, sd_intron, q_25, q_50, 
        q_75, q_95, q_99, q_999, q_9999, q_99999
    
    Functions:
        intron_stats:
            works with a dictionary returned by get_id_and_genes and compute introns functions
            makes a list of all intron lenghts from the genome intron set (acounts for repeated introns in different transcripts)
            uses that list to return numerical statistics about introns of a genome in a list with this order:
                [max_intron, min_intron, mean_intron,median_intron, sd_intron, q_25, q_50, 
                q_75, q_95, q_99, q_999, q_9999, q_99999]
        write_to_table: 
            uses taxonomy and extract modules as well as intron stats to get all the information 
            mentionaed in script summery for a gtf file
        
        A multiprocessor is used to run the write to table with 5 walkers and genomic_directory.csv as the input list

test_part_2.py (OUTDATED, I use notebook to make my graphes now!):
    Summery: uses the text file created by test_part_1.py to make multiple graphes with the information at hand
    
    makes the following figures:
    
    scatter with kbp and log10 bp scale of genome size vs max intron color coded for different kingdoms
    pie plots for distribution of kingdoms, assembly status, and assembly type of refseq data
    box plots for lenght distribution of introns in different quantiles across all genomes
        with kbp and log10 base pair scales
        with and without outliers

test_part_3.py 
        
    Summery: a script to find all the introns in our data by going through all the gtf files whose location is stored in a csv directory
    Multiprocessor is used to loop through all genomic files using get_id_and_introns function and record the output in introns.txt file

    
    Functions:
        get_id_and_introns: 
            extracts genome_id and genes from a genomic file, and returns the id and all the introns of the genome in a list.

test_part_4.py
    
    Summery: goes through all gtf files (gets their location from a csv directory) and exctracts information about their genes and gene introns

    Information it stores about each gene include: 
        gene_id genome_id kingdom phylum species gene_len non_coding_len ratio max_intron 
        max_intron_start max_intron_end total_sequence_length"
    
    The results is outputted in 4 different text files (because of using 4 workers):
        gene_1.txt gene_2.txt gene_3.txt gene_4.txt
    
    No analysis have been done on the table, the script is also still under furthur work ....

        
    
    
        
                
            
 
 