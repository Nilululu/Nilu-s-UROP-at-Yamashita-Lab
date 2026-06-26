Created on Thu Jun 25 16:45:26 2026


Hi, in this README I will give a overall review of my code structure, all the modules, 
and the order they work. I will also indicate if I need help or feedback on anything specific


    
Scripts I would like to get your review on:
    test_part_1.py 
        for this one specifically do you think I still struggle to understand multiprocessor 
        and put it to good use? if so what is it you think I haven't understood yet
    test_part_2.py


Parts I need help with:

    Some introns have a below zero length, if you find extract module, compute introns function, 
    you will see that I sort exons of each transcript by their start point before computing 
    the introns, I don't know why this happens unless the end of one exon is after the start
    of anotehr exon within the same transcript, is there a bug in my code, or does this have a 
    biological meaning?
    
    When I make box plots in test_part_2.py with log10 sacle, some of the intial boxplots 
    disapear (min, q_25, and q_50). You can see it in the figures saved in cluster, I don't know 
    what is causing it. I know there is no log10 for negetive numbers (in the case of intron lenghts calculated negetive), 
    but I am not sure if for example 50th quantile has negetive introns.

Scripst I recently run in the cluster:

    I run test_part_1.py and test_part_2.py
        you should be able to see the info_test2_logger.txt as well as all the figure save in png format 
        in cluster
    Some thoughts:
        for making figures, log10 or maybe any log seem to be much more informitive in graphes then kbp

Scripts I am trying to finish and run by the time we meet on Monday
    I want to make a test_part_3.py to make a violin plot of all the introns we have in our data 
    next steps would be to make the same plot for different kingdoms maybe to see if the distribution is 
    different between them
    
Goal out of this review and next week meetings:
    1. I want to know if the way I get data and make plots is correct and effcient and then move to asking more 
        questions
    2. Rob helped me create a document for questions that might be good to ask, 
        I know you have also added comments there. We can maybe talk about it 
    3. talking about the way we want to define a giant intron

ncbi_api.py:
    
    Summery: an script to download all reseq annotated genome folders from NCBI using a tsv file with 
    the accesion ids of genome folders and a API link

    files and folders created by running it: 
        ncbi_data_directory_0.1, genomic_directory.csv, error_dowload_log.txt, id_set.txt
    
    unzips and organize all successfully downloaded folders in ncbi_data_directory_0.1/
    
    logs warnings and errors in error_download_log.txt
    
    finds the location of the gtf file in the downlaoded folder and saves it in genomic_directory.csv
    
    keeps track of accession ids of all successfully downloaded genomic folders and save them in id_set.txt
    
    
extract.py:
    
    Summery: a module to take a gtf file as input and extract the information of interest about all genes in a genome 
    uses a dictionary data structure
    
    Data Structure:
        {gene_id: {"position": (start, end)", "transcrips":{}, "introns":set(), "exons":set(), "strand": "}}
        
            the "transcript" key has the follwing structure to its associated value
            {transcript_id: {"exons":[(start, end), ...] , "introns":[(start, end), ...]}}
    
    It works in two steps!
    To use it you should call extract_id_and_genes(gtf_file), this will give a dictionary with the data structure mentioned above
    Then call compute_introns on the said dictioanry, this will use the information about exons to compute and store information about introns
        
        
taxonomy.py:
    a module to take a taaxonomy number and use taxonomy dataset to return a taxonomy dictionary
    It needs names.dmp and nodes.dmp to be in th working directory 
    
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

test_part_2.py:
    Summery: uses the text file created by test_part_1.py to make multiple graphes with the information at hand
    
    makes the following figures:
    
    scatter with kbp and log10 bp scale of genome size vs max intron color coded for different kingdoms
    pie plots for distribution of kingdoms, assembly status, and assembly type of refseq data
    box plots for lenght distribution of introns in different quantiles across all genomes
        with kbp and log10 base pair scales
        with and without outliers
        
        
    
    
        
                
            
 
 