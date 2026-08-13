# BRCA2-999del5_cfDNA_nanopore

# Collection of scripts utilized as part of the (insert publication name here) study. 

## For the most part these scripts are wrappers and/or set up scripts utilized to prepare the Nanopore output for specific tools. Besides those there are several scripts primarily intended for statistical analysis and a few custom scripts for data extraction. All outside tools or scripts utilized for this analysis will be specified on this page.

The dependencies of each script will be detailed in their own dedicated directories.

The scripts that utilized BAM files were designed for files processed via Dorado which the official basecaller of Oxford Nanopore Technologies. Similarly scripts that utilize BED files were designed for files processed with Modkit pileup. 

Any input file types that are in a different format such as .csv and .tsv are derived from these two file types and the script utilized to create these files.

Example Dorado command line: 

time dorado-1.3.1-linux-x64/bin/dorado basecaller sup@v5.0.0 --kit-name SQK-NBD114-24 $tdir/pod5-cfDNA13-24 --modified-bases 5mCG_5hmCG -r --no-trim -v --device cuda:all > $tdir/myAoutputfiles13-24/FC13-24Goutput.bam
 
time dorado-1.3.1-linux-x64/bin/dorado demux --kit-name SQK-NBD114-24 $tdir/myAoutputfiles13-24/FC13-24Goutput.bam --output-dir $tdir/outputdemux13-24 --no-trim --emit-summary -v

time dorado-1.3.1-linux-x64/bin/dorado aligner hg38.fa $tdir/outputdemux13-24 -r --output-dir $tdir/myAoutputfiles13-24/output13-24 -v

Example Modkit command line: 
  time "$MODKIT_BIN" pileup "$bamfile" "$OUTPUT_DIR/${base_name}.bed" \
            --ref "$REF_GENOME" \
            --modified-bases 5mC 5hmC \
            --cpg \
            --combine-strands \
            --threads 8
        
