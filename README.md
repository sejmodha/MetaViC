# metaViC 
<html>
<head>
<h3>
metaViC : Virus metagenomics pipeline for unknown host or in absence of a host genome
</h3>
</head>

<body>
<p> We have developed an in-house metagenomics pipeline that we apply to the samples when we have sample from unknown hosts or if the host genome sequence is not sequenced. In the absence of the reference genome it is difficult to determine how to remove the host sequences from the sample. The pipeline described here can be applied directly to any virus metagenomics study. </p>
<p></p>
<p> This metagenomics pipeline is mainly divided into two major components. The first component deals with cleaning and removing non-viral contents from the reads and the second with de-novo assembly. </p>
<ol type="1">
  <li>Cleaning – removing known higher level organism and bacterial reads</li>
</ol>
<img align="left" src="images/Metagenomicspipeline1.png"/>
<p>
The first step of the cleaning pipeline is to remove the adapters using trim_galore (Krueger). In some cases of metagenomics sample preparation for sequencing ribosomal RNA (rRNA) removal might not be removed before sequencing the sample. In such cases and also in cases it is advisable to remove the rRNAs in silico using RiboPicker (Schmieder et al., 2012). After this, DIAMOND (Buchfink et al., 2014) can be run for each read against the refseq protein database for each file. Krona charts (Ondov et al., 2011) are generated for each DIAMOND output file that describes the read based classification of the sample. The DIAMOND results are then converted to a BLAST tabular output and the Genbank Identifier (GI) column from the output is extracted. GIs are mapped back to NCBI taxonomy databases to extract the corresponding taxonomy and division ID. In our case, we are only interested in the sequences that match viruses, environmental sequences and the sequences that do not match anything in the database. Therefore, any read that is matching protein sequences originating from bacteria, invertebrates, mammals, rodents, phages, plants, vertebrates, primates and synthetic constructs are identified based on the division ID and are filtered from the sample files. Finally, the sample files are properly paired using Prinseq-lite.pl (Schmieder et al., 2011). These reads can then be submitted to the next stage of the pipeline.
</p>

</body>
</html>

</body>
</html>

