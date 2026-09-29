<h1>Binary Locus Universal Extraction: Visualizing Fusion Genes</h1>
This application takes a BAM file submitted by a user, parses it and displays all fusion genes located. Fusion genes are displayed and color coded with a summary and additional information on each gene. Additionally, a viewer for chromosomes has been developed to allow the user to load in a human chromosome and display it in the Chrom-View. Both the fusion and chromosome viewers are able to be dragged left and right to view the rest of the reads along with their positions. More information on using the application along with the file structure is found below.<br></br>

* BLUE comes with sample BAM files to test with in the "Sample BAM" folder
* Application folder contains the following:
  * Python Applciation Scripts
    * bluefusion.py is the main script that generates the GUI and communicates with all subsequent modules
    * samfile.py parses a BAM file and returns relevant information such as fusion reads, contigs, etc.
  * Chromosomes folder
    * Place for the Homo_sapiens.GRCh38.dna.chromosome fils
    * They are too large to share and need to be downloaded and placed in this folder manually
    * Downloads
    * wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
    * wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{X,Y}.fa.gz
  * Genomes (GTF Reference file(s))
    * Place for the Homo_sapiens.GRCh38.100.db and Homo_sapiens.GRCh38.100.gtf
    * They are too large to share and need to be downloaded and placed in this folder manually
      * Download and unzip
      * wget ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
      * gunzip Homo_sapiens.GRCh38.100.gtf.gz

<h2>USING THE APPLICATION</h2>
<strong>Please note that this application leverages the 'pysam' library which no longer functions on Windows.</strong><br>
<strong>***Tested on Linux***</strong>

<h3>BAM ONLY:</h3>
To use this application, simply select a BAM file from the 'File' dropdown and the program will read it. If a fusion is found, the reads will be displayed in the Fusion Area in the middle of the application. The first fusion gene will be displayed and color coded on the left side, while the second will show up on the right. Positions mark the start of the base that follows and each 'marker' is 5 bases in length.

<h3>BAM W/ ANNOTATIONS:</h3>
In order for the program to retrieve the gene names and other relevant information, a GTF reference file will need to be loaded into the program. This can be done by selecting the reference the 'Ref' tab and the reference name associated with that genome. The program comes preloaded with Ensembles GChr38 annotation file for use. You can load the reference file before loading the BAM file and vice versa. For fusions found, the gene list A and the gene list B will populate with discovered fusions genes. More information about each gene can be found by clicking 'More Information' under each gene name. The summary section will show key information about the fusion that was found.

<h3>USING THE CHROM-VIEW:</h3>
The program allows each chromosome sequence to be viewed on a viewer under the Fusion-View. To load a chromosome sequence, select any of the chromosomes from the dropdown near the top of the window. The area below the Chromosome viewer updates in real-time the genes currently being viewed in the viewer and their positions. For the genes to be updated, the reference file must be loaded into the program, otherwise no genes will populate. Genes and positions can also be searched by using the search function next the the chromosome dropdown list. Type in the gene name or single position and the chrom-view will automatically go to the specified region/gene if it exists. The Fusion-View will update the same way, though will not show any sequences outside of the fusion sequences.

<h2>FUTURE IMPROVEMENTS</h2>

* QC/QA with other testing bam files - https://drive.google.com/drive/folders/1heeZSw72QU6n7sLH0qttaDlhdT4ePEkK?usp=sharing
  * Set 1 (BTF-474) original source - https://github.com/Oshlack/JAFFA/wiki/Example
  * Set 2 (genefuse) original source - http://opengene.org/dataset.html
* Add zoom feature to scale in and out of each base pair
* Add sequencing reads coverage plot for multi-reads view
* Show reads alignment direction in different color
* Add a summary table and allow users to download the report and graph
* Add links to external database such as FusionGDB and COSMIC
