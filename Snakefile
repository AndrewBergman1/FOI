import glob
import os

# Get the list of files in data/toUse/ (including extensions)
to_use_files = glob.glob("data/toUse/*")
original_files = [os.path.basename(f) for f in to_use_files]  # Includes extensions

# Get the list of datasets from data/cleaned_data/*.fsa
cleaned_files = glob.glob("data/cleaned_data/*.fsa")
datasets = [os.path.splitext(os.path.basename(f))[0] for f in cleaned_files]

rule all:
    input:
        # Cleaned databases
        expand("data/cleaned_data/{file}", file=original_files),

        # Concatenated database
        "data/cleaned_data/full_database.fasta",

        # Validated database
        "data/cleaned_data/validated_full_database.fasta",

        # CATCH probes
        #"outputs/catch/catch_probes.fna",

        # Syotti probes
        "outputs/syotti/syotti_0-baits.fna",
        "outputs/syotti/syotti_10-baits.fna",
        "outputs/syotti/syotti_20-baits.fna",
        "outputs/syotti/syotti_40-baits.fna",

        # Syotti kmer-counting files
        "outputs/kmer_freq/syotti/syotti.kmc_pre",
        "outputs/kmer_freq/syotti/syotti.kmc_suf",
        "outputs/kmer_freq/syotti/syotti_kmers.txt",

        # CATCH kmer-counting files
        #"outputs/kmer_freq/catch/catch.kmc_pre",
        #"outputs/kmer_freq/catch/catch.kmc_suf",
        #"outputs/kmer_freq/catch/catch_kmers.txt",

        # Dataset kmer-counting files
        expand("outputs/kmer_freq/datasets/{dataset}_kmers.txt", dataset=datasets),
        expand("outputs/kmer_freq/datasets/{dataset}.kmc_pre", dataset=datasets),
        expand("outputs/kmer_freq/datasets/{dataset}.kmc_suf", dataset=datasets),

        # Syotti UpSet plot
        # "/home/abergm/bait_gen/outputs/plots/syotti_upset_plot.png",

        # CATCH UpSet plot
        #"/home/abergm/bait_gen/outputs/plots/catch_upset_plot.png",

        # GC% and Tm plots
        "outputs/plots/0_gc_tm_boxplot.png",
        "outputs/plots/10_gc_tm_boxplot.png",
        "outputs/plots/20_gc_tm_boxplot.png",
        "outputs/plots/40_gc_tm_boxplot.png",

        #"outputs/plots/probe_sequence_coverage.png"


        # Kalamari database
        #"data/kalamari/kalamari.dbtype",

        # Indexed CATCH probes
        #"outputs/catch/catchDB",

        # Indexed Syotti probes
        "outputs/syotti/syotti_0_DB",
        "outputs/syotti/syotti_10_DB",
        "outputs/syotti/syotti_20_DB",
        "outputs/syotti/syotti_40_DB",

        
        # Indexed catch database
        "outputs/syotti/syotti_0_DB.idx",
        "outputs/syotti/syotti_10_DB.idx",
        "outputs/syotti/syotti_20_DB.idx",
        "outputs/syotti/syotti_40_DB.idx",

        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/syotti/syotti_alignments/0_mismatches/resultDB.1",
        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/syotti/syotti_alignments/10_mismatches/resultDB.1",
        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/syotti/syotti_alignments/20_mismatches/resultDB.1",
        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/syotti/syotti_alignments/40_mismatches/resultDB.1"
       # "/home/abergm/bait_gen/data/SILVA",
        "data/kalamari/kalamari.dbtype",
        "outputs/syotti/mtDNA_cpDNA_alignment/0/resultDB.1",
        "outputs/syotti/mtDNA_cpDNA_alignment/0/resultDB.m8",
        "outputs/syotti/mtDNA_cpDNA_alignment/0/problematic_probes.csv",
        "outputs/syotti/filtered_syotti_0-baits.fna",
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/human_genome.zip",
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/ncbi_dataset",
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/uppercase_GCF_000001405.26_GRCh38_genomic.fna",
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/human_genome_DB",
        "outputs/syotti/human_genome_alignment/0/resultDB.1",
        "outputs/syotti/human_genome_alignment/0/resultDB.m8",
        "outputs/syotti/human_genome_alignment/0/problematic_probes.csv",
        "outputs/syotti/complete_filtered_syotti_0-baits.fna",
        "outputs/plots/0_GC_Tm_filter_plots.png",
        "outputs/syotti/COMPLETED_0_SYOTTI_PROBES.fasta",
        "outputs/syotti/self_complementary/RNAfold_output.fold",
        "outputs/syotti/self_complementary/problematic_probes.csv",
        "/home/abergm/bait_gen/outputs/syotti/self_complementary/final_probes",
        "/home/abergm/bait_gen/outputs/plots/filter_comparison.png",
        "outputs/plots/free_energy_histogram.png",
        "outputs/kmer_freq/syotti/final/syotti.kmc_pre",
        "outputs/kmer_freq/syotti/final/syotti.kmc_suf",
        "outputs/kmer_freq/syotti/final/syotti_kmers.txt",
        "/home/abergm/bait_gen/outputs/plots/syotti_upset_plot.png",
        "outputs/plots/barchart_overlap.png"


        
        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/mitochondrial_search/mtDNA_cpDNA.fasta",
        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/mitochondrial_search/filtered_mtDNA_cpDNA.fasta",
        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/mitochondrial_search/db/mtDNA_cpDNA_DB",
        
        # "outputs/syotti/mtDNA_cpDNA_alignment/resultDB.1",
        # "outputs/syotti/mtDNA_cpDNA_alignment/resultDB.m8",
        # "outputs/syotti/mtDNA_cpDNA_alignment/problematic_probes.csv",
        # "outputs/syotti/filtered_syotti_10-baits.fna",
        
        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/human_genome.zip",
        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/ncbi_dataset",
        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/ncbi_dataset/data/GCF_000001405.26/uppercase_GCF_000001405.26_GRCh38_genomic.fna",
        # "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/ncbi_dataset/data/GCF_000001405.26/db/human_genome_DB",
        # "outputs/syotti/human_genome_alignment/resultDB.1",
        # "outputs/syotti/human_genome_alignment/resultDB.m8",
        # "outputs/syotti/human_genome_alignment/problematic_probes.csv",
        
        # "outputs/syotti/complete_filtered_syotti_10-baits.fna",
        # "outputs/plots/10_GC_Tm_filter_plots.png"  



        # Indexed syotti database
        #"outputs/syotti/syottiDB.idx",
        
        # Searching the NT database
        #"outputs/syotti/NT_alignment/resultDB.1",


        # CATCH probes alignment to kalamari database
        #"outputs/catch/kalamari_alignment/resultDB.1",

        # Syotti probes alignment to kalamari database
        #"outputs/syotti/kalamari_alignment/resultDB.1",
        
        # CATCH alignment in BLAST format
        #"outputs/catch/kalamari_alignment/resultDB.m8",


        # NT database
        #"/mnt/autofs/shared/storage01/users/exjobb/abergm/data/NT/ntDB"


#Clean datasets before merging
rule clean_datafiles:
    input:
        "data/toUse/{file}"
    output:
        "data/cleaned_data/{file}"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output})

        # Remove any blank lines from the input file and process with awk
        grep -v '^$' {input} | \
        awk '
        BEGIN {{
            RS=">"; ORS=""
        }}
        NR > 1 {{
            # Split the record into header and sequence
            n = split($0, lines, "\\n")
            header = lines[1]
            seq = ""
            for (i = 2; i <= n; i++) {{
                seq = seq lines[i]
            }}
            # Convert sequence to uppercase and remove non-ATCG characters
            seq = toupper(seq)
            gsub(/[^ATCG]/, "", seq)
            # Check sequence length and duplication
            if (length(seq) >= 120 && !seen[seq]++) {{
                print ">" header "\\n" seq "\\n"
            }}
        }}
        ' > {output}

        # Check if output file was created and is not empty
        if [ ! -s {output} ]; then
            echo "Error: Output file {output} was not created or is empty."
            exit 1
        fi
        """
        
# Concatenate the datasets to one large database
rule concat_datafiles:
    input:
        expand("data/cleaned_data/{file}", file=original_files)
    output:
        "data/cleaned_data/full_database.fasta"
    shell:
        """
        cat {input} > {output}
        """

# Ensure that the concatenated database doesnt contain any blank spaces or lowercase nucleotides
rule validate_database:
    input:
        "data/cleaned_data/full_database.fasta"
    output:
        "data/cleaned_data/validated_full_database.fasta"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output})

        # Remove any blank lines and validate sequences
        grep -v '^$' {input} | \
        awk '
        BEGIN {{
            RS = ">"
            ORS = ""
        }}
        NR > 1 {{
            # Remove leading and trailing whitespace
            gsub(/^\\s+|\\s+$/, "", $0)
            if ($0 != "") {{
                # Extract header and sequence
                n = split($0, lines, "\\n")
                header = lines[1]
                seq = ""
                for (i = 2; i <= n; i++) {{
                    seq = seq lines[i]
                }}
                # Remove whitespace and non-ATCG characters
                gsub(/[^ATCGatcg]/, "", seq)
                seq = toupper(seq)
                # Check if sequence length is at least 120
                if (length(seq) >= 120) {{
                    print ">" header "\\n" seq "\\n"
                }}
            }}
        }}
        ' > {output}

        # Check if output file was created and is not empty
        if [ ! -s {output} ]; then
            echo "Error: Output file {output} was not created or is empty."
            exit 1
        fi
        """

# Run syotti to generate probes
rule run_syotti:
    conda:
        "Syotti"
    input:
        "data/cleaned_data/validated_full_database.fasta"
    output:
        "outputs/syotti/syotti_{hd}-baits.fna",
        "outputs/syotti/syotti_{hd}-cover-fractions.txt",
        "outputs/syotti/syotti_{hd}-cover-marks.txt"
    shell:
        """
        bin/syotti/bin/syotti design --bait-len 120 --hamming-distance {wildcards.hd} -s {input} -r -o outputs/syotti/syotti_{wildcards.hd}
        """

# Create GC% and Tm boxplots for syotti probes
rule plot_gc_tm_boxplots:
    conda:
        "gc_tm_boxplot"
    input:
        #catch = "outputs/catch/catch_probes.fna",
        syotti = "outputs/syotti/syotti_{hd}-baits.fna"
    output:
        "outputs/plots/{hd}_gc_tm_boxplot.png"
    shell:
        "Rscript bin/gc_tm_boxplot.R {input.syotti} {output}"    





# Create syotti database to run against other databases
rule create_syotti_db:
    conda:
        "mmseqs2"
    input:
        "outputs/syotti/syotti_{hd}-baits.fna"
    output:
        "outputs/syotti/syotti_{hd}_DB"
    shell:
        "mmseqs createdb {input} {output}"

# Index the syotti databsae to run against other databsases
rule index_syotti_db: 
    conda:
        "mmseqs2"
    input:
        "outputs/syotti/syotti_{hd}_DB"
    output:
        "outputs/syotti/syotti_{hd}_DB.idx"
    shell:
        "mmseqs createindex {input} outputs/syotti/tmp --search-type 3"


# Download and index the SILVA database
rule download_index_SILVA: 
    conda:
        "mmseqs2"
    output:
        "/home/abergm/bait_gen/data/SILVA"
    params:
        tmpdir="/home/abergm/bait_gen/data/SILVA/tmp"
    shell:
        """
        mmseqs databases SILVA data/SILVA/SILVA {params.tmpdir}
        """

# Download and index the kalamari database
rule download_kalamari_database:
    conda:
       "mmseqs2"
    output:
       "data/kalamari/kalamari.dbtype"
    shell: 
       "mmseqs databases Kalamari data/kalamari/kalamari tmp"





# Download entries for chloroplasts and mitochondrial sequences to a fasta database 
rule search_mt:  
    conda:
        "entrez-direct"
    output: 
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/mitochondrial_search/mtDNA_cpDNA.fasta"
    shell:
        """
        esearch -db nucleotide -query "(mitochondrial[Title] OR chloroplast[Title]) AND (bacteria[Filter] OR fungi[Filter] OR protozoa[Filter])" | efetch -format fasta > {output} 
        """

# Filter the mtDNA/cpDNA database so that the entries are at least 120bp long
rule filter_mtDNA_cpDNA_database:
    input:
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/mitochondrial_search/mtDNA_cpDNA.fasta"
    output:
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/mitochondrial_search/filtered_mtDNA_cpDNA.fasta"
    shell:
        """
        bash scripts/filter_mtDNA_cpDNA.sh {input} {output}
        """


# Index the mtDNA/cpDNA database
rule index_mtDNA_cpDNA_database:
    conda:
        "mmseqs2"
    input: 
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/mitochondrial_search/filtered_mtDNA_cpDNA.fasta"
    output: 
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/mitochondrial_search/db/mtDNA_cpDNA_DB"
    shell:
        "mmseqs createdb {input} {output}"


# Align syotti probes (10 mismatches) to the cpDNA/mtDNA database
rule search_mtDNA_cpDNA_database:
    conda:
        "mmseqs2"
    input:
        syotti_db = "outputs/syotti/syotti_0_DB",
        mtDNA_cpDNA_db = "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/mitochondrial_search/db/mtDNA_cpDNA_DB"
    output:
        "outputs/syotti/mtDNA_cpDNA_alignment/0/resultDB.1"
    params:
        result_path = "outputs/syotti/mtDNA_cpDNA_alignment/0/resultDB",
        tmp_path = "outputs/syotti/mtDNA_cpDNA_alignment/0/tmp"
    shell:
        "mmseqs search {input.syotti_db} {input.mtDNA_cpDNA_db} {params.result_path} {params.tmp_path} --search-type 3 -s 7.0"
        


# Convert the alignment to blast format and make it easy to parse
rule convert_mtDNA_cpDNA_alignment_to_blast_format:
    conda:
        "mmseqs2"
    input:
        queryDB="outputs/syotti/syotti_0_DB",
        targetDB="/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/mitochondrial_search/db/mtDNA_cpDNA_DB",
        resultDB="outputs/syotti/mtDNA_cpDNA_alignment/0/resultDB.1"
    output:
        "outputs/syotti/mtDNA_cpDNA_alignment/0/resultDB.m8"
    params:
        db_prefix="outputs/syotti/mtDNA_cpDNA_alignment/0/resultDB"
    shell:
        "mmseqs convertalis {input.queryDB} {input.targetDB} {params.db_prefix} {output} --format-output 'query,target,bits,pident,fident,nident'"


# Parse the blast-format alignment file, identify probes that align to mtDNA/cpDNA with identity > 0.5
rule identify_problematic_probes:
    input:
        "outputs/syotti/mtDNA_cpDNA_alignment/0/resultDB.m8"
    output:
        "outputs/syotti/mtDNA_cpDNA_alignment/0/problematic_probes.csv"
    shell:
        "python3 scripts/filter_problematic_probes.py {input} {output}"


# Remove the probes that align with >0.5 identity
rule remove_problematic_probes:
    input:
        baits = "outputs/syotti/syotti_0-baits.fna",
        problematic_baits = "outputs/syotti/mtDNA_cpDNA_alignment/0/problematic_probes.csv"
    output:
        "outputs/syotti/filtered_syotti_0-baits.fna"
    shell:
        "python3 scripts/remove_problematic_probes.py {input.baits} {input.problematic_baits} {output}"


# Download the human genome (GRCh38)
rule download_human_genome: 
    output:
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/human_genome.zip"
    shell:
        """
        wget -O /mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/human_genome.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.26/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
        """

# Unzip the human genome 
rule unzip_human_genome:
    input: 
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/human_genome.zip"
    output:
        directory("/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/ncbi_dataset")
    shell:
        "unzip {input}"


# Ensure uppercase database
rule process_human_genome: 
    input: 
        req = "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/ncbi_dataset",
        hum_genome = "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/ncbi_dataset/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna"
    output: 
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/uppercase_GCF_000001405.26_GRCh38_genomic.fna"
    shell: 
        """ 
        cat {input.hum_genome} | tr "[atcg]" "[ATCG]" > {output}
        """

# Index the human genome for alignment
rule index_human_genome: 
    conda:
        "mmseqs2"
    input:
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/uppercase_GCF_000001405.26_GRCh38_genomic.fna"
    output:
        "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/human_genome_DB"
    shell:
        "mmseqs createdb {input} {output}"


# Align syotti probes to the human genome
rule search_human_genome:
    conda:
        "mmseqs2"
    input:
        syotti_db = "outputs/syotti/syotti_0_DB",
        human_genome = "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/human_genome_DB"

    output:
        "outputs/syotti/human_genome_alignment/0/resultDB.1"
    params:
        result_path = "outputs/syotti/human_genome_alignment/0/resultDB",
        tmp_path = "outputs/syotti/human_genome_alignment/0/tmp"
    shell:
        "mmseqs search {input.syotti_db} {input.human_genome} {params.result_path} {params.tmp_path} --search-type 3 -s 7.0"


# Convert alignment to blast format for easy parsing
rule convert_human_alignment_to_blast_format:
    conda:
        "mmseqs2"
    input:
       queryDB = "outputs/syotti/syotti_0_DB",
       targetDB = "/mnt/autofs/shared/storage01/users/exjobb/abergm/bait_gen/human_genome_search/human_genome_DB",
       #resultDB = "outputs/syotti/kalamari_alignment/0/resultDB.1" 
    output:
       "outputs/syotti/human_genome_alignment/0/resultDB.m8"
    params:
        db_prefix="outputs/syotti/human_genome_alignment/0/resultDB"
    shell:
        """
        mmseqs convertalis {input.queryDB} {input.targetDB} {params.db_prefix} {output} --format-output "query,target,bits,pident,fident,nident"
        """   

# Identify probes that align with >0.8 identity
rule identify_human_problematic_probes:
    input:
        "outputs/syotti/human_genome_alignment/0/resultDB.m8"
    output:
        "outputs/syotti/human_genome_alignment/0/problematic_probes.csv"
    shell:
        "python3 scripts/filter_problematic_probes.py {input} {output}"

# Remove probes that align to human genome with >0.8 identity
rule remove_human_problematic_probes: 
    input:
        baits = "outputs/syotti/filtered_syotti_0-baits.fna",
        problematic_baits = "outputs/syotti/human_genome_alignment/problematic_probes.csv"
    output:
        "outputs/syotti/complete_filtered_syotti_0-baits.fna"
    shell:
        "python3 scripts/remove_problematic_probes.py {input.baits} {input.problematic_baits} {output}"
       
rule filter_gc:
    conda:
        "filter_gc"
    input:
        baits = "outputs/syotti/complete_filtered_syotti_0-baits.fna"
    output:
        plot = "outputs/plots/0_GC_Tm_filter_plots.png",
        probes = "outputs/syotti/COMPLETED_0_SYOTTI_PROBES.fasta"
    shell:
        "Rscript bin/filter_gc.R {input.baits} 0"

rule find_self_complementarity: 
    conda:
        "viennarna"
    input:
        "cdCOMPLETED_0_SYOTTI_PROBES.fasta"
    output:
        "outputs/syotti/self_complementary/RNAfold_output.fold"
    params:
        "outputs/syotti/self_complementary"
    shell:
        """
        mkdir -p {params}
        cd {params}
        RNAfold --infile={input} --outfile=RNAfold_output.fold --jobs=16
        """

rule find_problematic_free_energy:
    conda:
        "filter_free_energy"
    input:
        "outputs/syotti/self_complementary/RNAfold_output.fold"
    output:
        probes = "outputs/syotti/self_complementary/problematic_probes.csv",
        plot = "outputs/plots/free_energy_histogram.png"
    shell:
        """
        Rscript bin/filter_free_energy.R {input} {output.probes} {output.plot}
        """

rule remove_problematic_free_energy:
    conda:
        "filter_gc" # This en has the required packages already
    input:
        problematic_probes = "outputs/syotti/self_complementary/problematic_probes.csv",
        probes = "/home/abergm/bait_gen/outputs/syotti/COMPLETED_0_SYOTTI_PROBES.fasta",
    output:
        "/home/abergm/bait_gen/outputs/syotti/self_complementary/final_probes"
    shell:
        """
        Rscript bin/remove_self_compl.R {input.probes} {input.problematic_probes} {output}
        """

rule plot_probe_diff:
    conda:
        "filter_gc"
    input:
        original_probes = "/home/abergm/bait_gen/outputs/syotti/syotti_0-baits.fna",
        filtered_probes = "/home/abergm/bait_gen/outputs/syotti/self_complementary/final_probes",
        semifiltered_probes = "/home/abergm/bait_gen/outputs/syotti/COMPLETED_0_SYOTTI_PROBES.fasta"
    output:
        "/home/abergm/bait_gen/outputs/plots/filter_comparison.png"
    shell:
        """
        Rscript bin/probe_filter_comparison.R {input.original_probes} {input.filtered_probes} {input.semifiltered_probes} {output}
        """

rule count_syotti_kmers:
    conda:
        "kmc"
    input:
        "outputs/syotti/self_complementary/final_probes"
    output:
        kmc_pre = "outputs/kmer_freq/syotti/final/syotti.kmc_pre",
        kmc_suf = "outputs/kmer_freq/syotti/final/syotti.kmc_suf",
        kmers_txt = "outputs/kmer_freq/syotti/final/syotti_kmers.txt"
    params:
        tmp_dir = "outputs/kmer_freq/syotti/final/tmp"
    shell:
        """
        mkdir -p {params.tmp_dir}
        kmc -k120 -fa -ci1 {input} outputs/kmer_freq/syotti/final/syotti {params.tmp_dir}
        kmc_tools transform outputs/kmer_freq/syotti/final/syotti dump {output.kmers_txt}
        """


rule plot_upset_syotti:
    conda:
        "kmer_overlap"
    input: 
        "outputs/kmer_freq/syotti/final/syotti_kmers.txt",    
    output:
        plot = "/home/abergm/bait_gen/outputs/plots/syotti_upset_plot.png",
        table = "outputs/tables/kmer_overlap_post_filtering.tsv"
    params:
        table_dir = "outputs/tables"
    shell:
        """
        mkdir -p {params.table_dir}
        Rscript bin/plot_syotti_kmer_overlap.R {output.table}
        """

rule generate_overlap_barchart:
    conda:
        "kmer_overlap"
    input:
        table = "outputs/tables/kmer_overlap_post_filtering.tsv"
    output:
        plot = "outputs/plots/barchart_overlap.png"
    shell:
        "Rscript bin/generate_overlap_plot.R {input.table} {output.plot}"





# rule search_syotti_hits:
#     conda:
#         "mmseqs2"
#     input:
#         syottiDB = "outputs/syotti/syotti_{hd}_DB",
#         NT = "/mnt/autofs/shared/storage01/users/exjobb/abergm/data/NT/ntDB"
#     params:
#         tmp = "/mnt/autofs/shared/storage01/users/exjobb/abergm/syotti/syotti_alignments/{hd}_mismatches/tmp"
#     output:
#         "/mnt/autofs/shared/storage01/users/exjobb/abergm/syotti/syotti_alignments/{hd}_mismatches/resultDB.1"
#     shell:
#         "mmseqs search {input.syottiDB} {input.NT} outputs/syotti/NT_alignment/{wildcards.hd}_mismatches/resultDB {params.tmp} --search-type 3 -s 1.0 -a"



# rule count_dataset_kmers:
#     conda:
#         "kmc"
#     input:
#         "data/cleaned_data/{dataset}.fsa"
#     output:
#         kmc_pre = "outputs/kmer_freq/datasets/{dataset}.kmc_pre",
#         kmc_suf = "outputs/kmer_freq/datasets/{dataset}.kmc_suf",
#         kmers_txt = "outputs/kmer_freq/datasets/{dataset}_kmers.txt"
#     params:
#         tmp_dir = "outputs/kmer_freq/datasets/{dataset}_tmp"
#     shell:
#         """
#         mkdir -p {params.tmp_dir}
#         kmc -k120 -fa -ci1 {input} outputs/kmer_freq/datasets/{wildcards.dataset} {params.tmp_dir}
#         kmc_tools transform outputs/kmer_freq/datasets/{wildcards.dataset} dump {output.kmers_txt}
#         """

# rule run_catch:
#     conda:
#         "catch"
#     input:
#         "data/cleaned_data/validated_full_database.fasta"
#     output:
#         "outputs/catch/catch_probes.fna"
#     shell:
#         """
#         python3 bin/catch/catch/bin/design.py {input} -o {output} --island-of-exact-match 120 -pl 120 --verbose
#         """


# rule count_catch_kmers:
#     conda:
#         "kmc"
#     input:
#         "outputs/catch/catch_probes.fna"
#     output:
#         kmc_pre = "outputs/kmer_freq/catch/catch.kmc_pre",
#         kmc_suf = "outputs/kmer_freq/catch/catch.kmc_suf",
#         kmers_txt = "outputs/kmer_freq/catch/catch_kmers.txt"
#     params:
#         tmp_dir = "outputs/kmer_freq/catch/tmp"
#     shell:
#         """
#         mkdir -p {params.tmp_dir}
#         kmc -k120 -fa -ci1 {input} outputs/kmer_freq/catch/catch {params.tmp_dir}
#         kmc_tools transform outputs/kmer_freq/catch/catch dump {output.kmers_txt}
#         """







# rule plot_sequence_coverage:
#     conda:
#         "plot_seq_coverage"
#     input:
#         probes = "outputs/syotti/syotti-baits.fna",
#         fasta = "data/cleaned_data/validated_full_database.fasta"
#     output: 
#         "outputs/plots/probe_sequence_coverage.png"
#     shell: 
#         "Rscript bin/plot_seq_coverage.R {input.fasta} {input.probes}"

# rule plot_upset_catch:
#     conda:
#         "kmer_overlap"
#     input:
#         "outputs/kmer_freq/catch/catch_kmers.txt",
#         "outputs/kmer_freq/syotti/syotti_kmers.txt"
#     output:
#         "/home/abergm/bait_gen/outputs/plots/catch_upset_plot.png"
#     shell:
#         "Rscript bin/plot_catch_kmer_overlap.R"

# rule create_catch_db:
#     conda:
#         "mmseqs2"
#     input:
#         "outputs/catch/catch_probes.fna"
#     output:
#         "outputs/catch/catchDB"
#     shell:
#         "mmseqs createdb {input} {output}"


# rule index_catch_db: 
#     conda:
#         "mmseqs2"
#     input:
#         "outputs/catch/catchDB"
#     output:
#         "outputs/catch/catchDB.idx"
#     shell:
#         "mmseqs createindex {input} outputs/catch/tmp --search-type 3"




#rule catch_search_kalamari_taxa: 
#    conda:
#        "mmseqs2"
#    input:
#        catchDB = "outputs/catch/catchDB",
#        kalamari = "data/kalamari/kalamari"
#    output:
#        "outputs/catch/kalamari_alignment/resultDB.1"
#    shell:
#        "mmseqs search {input.catchDB} {input.kalamari} outputs/catch/kalamari_alignment/resultDB outputs/catch/kalamari_alignment/tmp --search-type 3 -s 7.0"


#rule convert_catch_alignment_to_blast_format:
#    conda:
#        "mmseqs2"
#    input:
#        queryDB = "outputs/catch/catchDB",
#        targetDB = "data/kalamari/kalamari",
#        resultDB = "outputs/catch/kalamari_alignment/resultDB.1"
#    output: 
#        "outputs/catch/kalamari_alignment/resultDB.m8"
#    shell:
#        """
#        mmseqs convertalis {input.queryDB} {input.targetDB} outputs/catch/kalamari_alignment/resultDB {output} --format-output "query,target,taxlineage,taxname,bits"
#        """




# rule plot_catch_taxa 
# rule plot_syotti_taxa 