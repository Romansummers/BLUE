import pysam
import tkinter as tk

from pyensembl import Genome, fasta

ALL_QUERIES = []
ALL_POSITIONS = []
ALL_LENGTHS = []
ALL_CHROMOSOMES = []
ALL_SEQUENCES = []
ALL_CIGARS = []

FUSION_POSITIONS = []
FUSION_LENGTHS = []


class InitFiles:
    def __init__(self, filepath, filetype):
        self.filepath = filepath
        self.filetype = filetype
        self.reads = []

    def init_chr_file(self):
        try:
            self.filetype = self.filetype.replace("CHR", "")  # To access the FASTA, we need the exact chromosome name
            ChromosomeLoad().load_chromosome(self.filepath, self.filetype)
        except Exception:
            raise Exception("Error loading chromosome")
        finally:
            root.quit()
            root.destroy()

    def init_sam_file(self):
        try:
            self.loading_msg("\n\nLoading BAM file...")
            BAMParser().load_sam(self.filepath)
            self.loading_msg("\n\nStoring BAM data...")
            BAMParser().store_all_reads()
            self.loading_msg("\n\nLocating fusions...")
            self.reads = BAMParser().map_fusion_reads()
        except Exception:
            raise Exception("Error loading bam file")
        finally:
            root.quit()
            root.destroy()

    def init_gtf_file(self):
        try:
            GenomeLoad().load_genome(self.filepath)
        except Exception:
            raise Exception("Error loading the gtf file")
        finally:
            root.quit()
            root.destroy()

    def create_loading(self):
        global root
        root = tk.Tk()
        root.title("Loading")
        root.geometry("250x100")
        try:
            if "CHR" in self.filetype:
                label = tk.Label(root, text="\n\nLoading chromosome...\n")
                label.pack()
                root.after(200, self.init_chr_file)
            if self.filetype == "SAM":
                label = tk.Label(root, text="\n\nLoading BAM file...\n")
                label.pack()
                root.after(200, self.init_sam_file)
            if self.filetype == "GTF":
                label = tk.Label(root, text="\n\nLoading/creating reference database..")
                label.pack()
                root.after(200, self.init_gtf_file)
        except Exception:
            raise Exception

        root.resizable(False, False)
        root.mainloop()

    @staticmethod
    def loading_msg(msg):
        var1 = tk.StringVar()
        label = tk.Label(root, textvariable=var1)
        label.pack()
        var1.set(msg)


class ChromosomeLoad:
    @staticmethod
    def load_chromosome(file,type):
        global chromdata
        chromdata = fasta.parse_fasta_dictionary(file)[type]


class GenomeLoad:
    # Loads reference genome and stores it in a glabal variable
    @staticmethod
    def load_genome(file):
        global genedata
        genedata = Genome(reference_name='GRCh38', annotation_name='ENSEMBL', gtf_path_or_url=file)
        genedata.index()


class BAMParser:
    # Loads file and stores it in a global variable for class to use
    def load_sam(self, file):
        global results

        sam = pysam.AlignmentFile(file, "rb")
        results = sam.fetch()

    # Stores every read in list to be referenced by index
    @staticmethod
    def store_all_reads():
        global ALL_QUERIES
        global ALL_POSITIONS
        global ALL_LENGTHS
        global ALL_CHROMOSOMES
        global ALL_SEQUENCES
        global ALL_CIGARS
        ALL_QUERIES = []
        ALL_POSITIONS = []
        ALL_LENGTHS = []
        ALL_CHROMOSOMES = []
        ALL_SEQUENCES = []
        ALL_CIGARS = []

        for read in results:
            # Format contig to be numeric
            chromosome = read.reference_name
            chromosome = chromosome.replace('chr', '')
            chromosome = int(chromosome)

            ALL_QUERIES.append(read.query_name)
            ALL_POSITIONS.append(read.pos + 1)  # pysam counts start position one behind the position in the bam file
            ALL_LENGTHS.append(read.query_length)
            ALL_CHROMOSOMES.append(chromosome)
            ALL_SEQUENCES.append(read.seq)
            ALL_CIGARS.append(read.cigar)

    # Using the data stored, locate any fusion reads within the bam file and show it on screen
    def map_fusion_reads(self):
        global FUSION_POSITIONS
        global FUSION_LENGTHS

        FUSION_READS = []
        FUSION_POSITIONS = []
        FUSION_LENGTHS = []

        positions = self.get_unique_positions()
        fusion_position = 0
        supp_position = 0
        position_headers = []  # Initialize positions for the reads
        paired_reads = []

        i = len(positions)
        right_indx, left_indx = 0, 0
        total_fusions = 0  # Stores the number of fusion reads found
        isFirst = False  # Used to create headers if fusion was found
        isSet = False  # Used to capture the first index for comparisons against all subsequent ones
        start = 0  # Index incrememented after each comparison
        x = 0

        # Initialize structure of dictionary that will store fusion read information
        read_dict = {"fusion": "", "fusionName": "", "fusionId": "",
                     "fusionPosition": fusion_position, "fusionLength": "",
                     "supplement": "", "supplementName": "", "supplementId": "",
                     "supplementPosition": supp_position, "supplementLength": "",
                     "isFirst": isFirst, "headers": position_headers}

        if i > 1:
            while x < len(ALL_POSITIONS):
                right_indx = x

                if not isSet:
                    left_indx = x
                    isSet = True

                x = x + 1
                # If we reach a new instance of positions, we want to test them to see if the reads fuse
                if ALL_POSITIONS[left_indx] != ALL_POSITIONS[right_indx] and left_indx not in paired_reads and right_indx not in paired_reads:
                    read_dict['fusion'], fusion_indx, read_dict['supplement'], supp_indx = self.get_fusion_read(left_indx, right_indx)
                    paired_reads.append(left_indx)
                    paired_reads.append(right_indx)

                    if read_dict['fusion'] != "No Fusions Detected":
                        total_fusions = total_fusions + 1
                        # Store the fusion reads for larger scope reference
                        FUSION_POSITIONS.append(ALL_POSITIONS[fusion_indx])
                        FUSION_LENGTHS.append(ALL_POSITIONS[fusion_indx])

                        read_dict['fusionLength'] = ALL_LENGTHS[fusion_indx]
                        read_dict['fusionPosition'] = ALL_POSITIONS[fusion_indx]
                        fusion_chr = ALL_CHROMOSOMES[fusion_indx]

                        read_dict['supplementLength'] = ALL_LENGTHS[supp_indx]
                        read_dict['supplementPosition'] = ALL_POSITIONS[supp_indx]
                        supp_chr = ALL_CHROMOSOMES[supp_indx]

                        if total_fusions == 1:
                            read_dict['isFirst'] = True
                            read_dict['headers'] = self.get_position_headers(read_dict['fusionLength'], read_dict['fusionPosition'], 50)

                        # If no reference genome has been loaded, we only can display the sequences, not names
                        try:
                            genedata
                        except NameError:
                            read_dict['fusionName'], read_dict['fusionId'], read_dict['supplementName'], read_dict['supplementId'] = "", "", "", ""
                        else:
                            read_dict['fusionName'] = self.get_gene_names(read_dict['fusionLength'], read_dict['fusionPosition'], fusion_chr)
                            read_dict['fusionId'] = self.get_gene_ids(read_dict['fusionLength'], read_dict['fusionPosition'], fusion_chr)[0]

                            read_dict['supplementName'] = self.get_gene_names(read_dict['supplementLength'], read_dict['supplementPosition'], supp_chr)
                            read_dict['supplementId'] = self.get_gene_ids(read_dict['supplementLength'], read_dict['supplementPosition'], supp_chr)[0]

                    if read_dict['fusion'] != "No Fusions Detected":
                        FUSION_READS.append(dict(read_dict))  # Only append valid fusion reads
                    start = start + 1
                    x = start
                    isSet = False
                    read_dict['isFirst'] = False

        if total_fusions == 0:
            FUSION_READS.append(dict(read_dict))  # If no fusions found display one instance of alert

        return FUSION_READS

    # Checks if two reads are a fusion based on unique indeces
    @staticmethod
    def get_fusion_read(left_indx, right_indx):
        fusion = "No Fusions Detected"
        fusion_indx = -1
        supplement = "No Fusions Detected"
        supp_indx = -1

        Read1Length = ALL_LENGTHS[left_indx]
        Read1Seq = ALL_SEQUENCES[left_indx]

        Read2Length = ALL_LENGTHS[right_indx]
        Read2Seq = ALL_SEQUENCES[right_indx]

        if (Read1Length > Read2Length) and (Read1Seq.endswith(Read2Seq) or Read1Seq.startswith(Read2Seq)):
            fusion = Read1Seq
            fusion_indx = left_indx
            supplement = Read2Seq
            supp_indx = right_indx
        elif (Read2Length > Read1Length) and (Read2Seq.endswith(Read1Seq) or Read2Seq.startswith(Read1Seq)):
            fusion = Read2Seq
            fusion_indx = right_indx
            supplement = Read1Seq
            supp_indx = left_indx

        return fusion, fusion_indx, supplement, supp_indx

    # This takes a number and returns an evenly distributed list of positions to use as headers
    @staticmethod
    def get_position_headers(length, start_position, pos_between):
        positions = []

        for x in range(length):
            if x == 0:
                positions.append(start_position)
            if x % pos_between == 0 and x != 0:  # This states the number of base pairs before a position header
                position = start_position + x
                positions.append(position)
            elif x == length-1:
                position = start_position + length
                positions.append(position)
        return positions

    # Use this to reference the found fusions
    @staticmethod
    def get_fusions():
        return FUSION_POSITIONS, FUSION_LENGTHS

    # Get the unique reads
    @staticmethod
    def get_unique_positions():
        unique_positions = []
        for pos in ALL_POSITIONS:
            if pos not in unique_positions:
                unique_positions.append(pos)

        return unique_positions

    # Get the gene names associates with each read
    @staticmethod
    def get_gene_names(length, position, chromosome):
        genes = []
        start_position = position
        end_position = start_position + length

        gene_names = genedata.genes_at_locus(contig=chromosome, position=start_position, end=end_position)
        for gene in gene_names:
            if gene.biotype == 'protein_coding':
                genes.append(gene.gene_name)

        return genes

    # Get the ID's of genes based on size, position and chromosome
    @staticmethod
    def get_gene_ids(length, position, chromosome):
        genes = []
        start_position = position
        end_position = start_position + length
        chromosome = chromosome

        gene_names = genedata.genes_at_locus(contig=chromosome, position=start_position, end=end_position)
        for gene in gene_names:
            if gene.biotype == 'protein_coding':
                genes.append(gene.gene_id)

        return genes

    # Gets all genes associated with a BAM file
    def get_all_gene_names(self):
        read_gene_names = []

        for i in range(len(ALL_POSITIONS)):
            read_gene_names[i] = self.get_gene_names(ALL_LENGTHS[i], ALL_POSITIONS[i], ALL_CHROMOSOMES[i])
        return read_gene_names

    # Given a position and chromosome, return info about all genes located
    @staticmethod
    def get_gene_info(start_position, chromosome):
        gene_names = genedata.genes_at_locus(contig=chromosome, position=start_position)
        genes = []
        for gene in gene_names:
            genes.append(gene)
        return genes

    # Given a gene name, return it's ID
    @staticmethod
    def get_id_by_name(gene_name):
        gene_id = genedata.gene_ids_of_gene_name(gene_name)
        return gene_id[0]

    # Returns a string position range of a gene given a name (with size)
    @staticmethod
    def get_positions_by_name(gene_name):
        positions = genedata.loci_of_gene_names(gene_name)
        start = positions[0].start
        end = positions[0].end
        return str(start), str(end)

    # Returns the starting position of a gene by name
    @staticmethod
    def get_start_position_by_name(gene_name):
        positions = genedata.loci_of_gene_names(gene_name)
        start_position = positions[0].start
        return start_position

    # Returns the chromosome of a given gene name
    @staticmethod
    def get_chromosome_by_name(gene_name):
        gene = genedata.genes_by_name(gene_name)
        chromsome = gene[0].contig
        return chromsome

    # Get the transcript ids using the gene name
    @staticmethod
    def get_transcript_ids_by_name(gene_name):
        transcript_ids = []
        transcripts = genedata.transcript_ids_of_gene_name(gene_name)
        for t in transcripts:
            transcript_ids.append(t)
        transcript_ids = sorted(transcript_ids)
        return transcript_ids

    # Get the name of the transcript using the transcript id
    @staticmethod
    def get_transcript_name_by_id(t_id):
        transcript = []
        t_name = genedata.transcript_name_of_transcript_id(t_id)
        transcript.append(t_name)
        return transcript[0]

    # Get the start and end position of a transcript given an id
    @staticmethod
    def get_transcript_positions_by_id(t_id):
        positions = genedata.locus_of_transcript_id(t_id)
        start = positions.start
        end = positions.end
        return str(start), str(end)

    # Get all exon regions associated with a transcript id
    @staticmethod
    def get_exon_regions_by_transcript_id(t_id):
        exon_starts = []
        exon_ends = []
        exon_ids = genedata.exon_ids_of_transcript_id(t_id)
        for exon_id in exon_ids:
            start = genedata.locus_of_exon_id(exon_id).start
            end = genedata.locus_of_exon_id(exon_id).end
            if start not in exon_starts:
                exon_starts.append(start)
                exon_ends.append(end)

        exon_starts = sorted(exon_starts)
        exon_ends = sorted(exon_ends)
        return exon_starts, exon_ends

    # Returns exon regions given a gene name
    @staticmethod
    def get_exon_regions_by_name(gene_name):
        exon_regions = []
        exon_ids = genedata.exon_ids_of_gene_name(gene_name)
        for exon_id in exon_ids:
            start = genedata.locus_of_exon_id(exon_id).start
            end = genedata.locus_of_exon_id(exon_id).end
            exon_string = str(start) + " - " + str(end)
            if exon_string not in exon_regions:
                exon_regions.append(exon_string)

        exon_regions = sorted(exon_regions)
        return exon_regions

    # Returns chromosome sequence
    @staticmethod
    def get_chromosome_seq():
        return chromdata

    # Return the position of a sequence in the loaded chromosome sequence
    @staticmethod
    def get_chrom_seq_pos(seq):
        return chromdata.find(seq)

    # Stores the cigar in two lists to be referenced by index
    # 0 = match; 1 = insertion; 2 = deletion; 3 = skip
    # 4 = soft clipping; 5 = hard clipping, 6 = padding
    @staticmethod
    def parse_cigar(cigar):
        cigar_type = []
        cigar_length = []
        for (cigarType, cigarLength) in cigar:
            cigar_type.append(cigarType)
            cigar_length.append(cigar_length)
        return cigar_type, cigar_length
