import tkinter as tk
import webbrowser

from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox

from samfile import InitFiles, BAMParser

#########################################
# Files required:
# Annotation file
# ****ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
# BAM file & Index file
#########################################

VIEWER_SPAN = 0  # Positional span of the chromosome viewer (e.g. 10,000 bases)
CHROM_START = 0  # Beginning position of the chromosome viewer
CHROM_END = 0  # End position of the chromosome viewer
CHROMOSOME_SEQ = ""  #
SEL_CHROMOSOME = ""
GENOME_LOADED = False
BAM_LOADED = False

FUSION_GENES = []

DEFAULT_THEME = "lightsteelblue1"  # The color of the application


class FusionGUI:
    """
    The Fusion GUI is initialized using the runblue method.
    This application is used in conjuction with samfile.py which
    parses bam files and displays any findings on the GUI.

    NOTE: samfile import uses pysam library to parse BAM files which
    has been depricated on Windows and will run on Linux and MacOS.

    All widgets are prone to consistent change and so they are
    be opened to all other classes and methods below this one.

    The widgets are stored in a list to be referenced by classes
    that deal with the GUI as a whole (e.g. color changes).

    Within the mainloop the application checks consistently for
    any genes to be displayed in the chromosome viewer.
    NOTE: Both a chromosome AND a reference genome need to be
    loaded in order to:
    a) Search by gene name
    b)Display gene names when using the chromosome viewer

    """

    def runblue(self):
        Search = SearchFunctions(self)
        MenuFunc = MenuFunctions(self)
        Mouse = MouseFunctions(self)

        self.win = Tk()
        self.win.title("BLUE")
        self.win.geometry("1000x750")  # Size of window
        self.win.config(bg=DEFAULT_THEME)

        # Create a menu bar
        menu_bar = Menu(self.win)
        self.win.config(menu=menu_bar)

        # Create menu and add menu items
        file_menu = Menu(menu_bar, tearoff=0)
        file_menu.add_command(label='Load from BAM file...', command=MenuFunc.load_bam_file)
        file_menu.add_separator()
        file_menu.add_command(label='Exit', command=self._quit)
        menu_bar.add_cascade(label='File', menu=file_menu)
        # Reference human genome menu
        hg_menu = Menu(menu_bar, tearoff=0)
        hg_menu.add_command(label='hg38', command=lambda: MenuFunc.load_genome_file("hg38"))
        menu_bar.add_cascade(label='Ref', menu=hg_menu)
        # Theme options
        theme_menu = Menu(menu_bar, tearoff=0)
        theme_menu.add_command(label="Default", command=lambda: MenuFunc.change_theme("lightsteelblue1"))
        theme_menu.add_command(label="Lavender", command=lambda: MenuFunc.change_theme("lavender"))
        theme_menu.add_command(label="Rose", command=lambda: MenuFunc.change_theme("misty rose"))
        theme_menu.add_command(label="Pinkesque", command=lambda: MenuFunc.change_theme("thistle2"))
        theme_menu.add_command(label="Antique", command=lambda: MenuFunc.change_theme("antique white"))
        theme_menu.add_command(label="Sky", command=lambda: MenuFunc.change_theme("lightskyblue1"))
        theme_menu.add_command(label="Aqua", command=lambda: MenuFunc.change_theme("PaleTurquoise2"))
        theme_menu.add_command(label="Snow White", command=lambda: MenuFunc.change_theme("snow2"))
        theme_menu.add_command(label="Graydation", command=lambda: MenuFunc.change_theme("gray81"))
        menu_bar.add_cascade(label="Theme", menu=theme_menu)

        help_menu = Menu(menu_bar, tearoff=0)
        help_menu.add_command(label="How to...", command=MenuFunc.open_how_to)
        help_menu.add_separator()
        help_menu.add_command(label="Report to admin...", command=MenuFunc.report_error)
        menu_bar.add_cascade(label="Help", menu=help_menu)

        self.frameRef = Frame(self.win)
        self.frameRef.pack(side=TOP)

        self.refLabel = Label(self.frameRef, text='No reference in Use', bg=DEFAULT_THEME)
        self.refLabel.grid(row=0, column=0)

        self.frame1 = Frame(self.win, bg=DEFAULT_THEME)
        self.frame1.pack()

        # Add a label; enter a gene or a locus
        self.searchlbl1 = Label(self.frame1, text='Enter a gene or position:', bg=DEFAULT_THEME)
        self.searchlbl1.grid(column=3, row=0)

        # Adding a Text box Entry widget for enter a gene or a locus
        search_var = tk.StringVar()
        search_box = Entry(self.frame1, width=12, textvariable=search_var)
        search_box.grid(column=4, row=0)

        # Add search button for a gene or a locus
        search_button = Button(self.frame1, text='Search', command=lambda: Search.search_position(search_box.get()))
        search_button.grid(column=5, row=0)

        # Add a label with Chromosome drop-down
        self.chrlbl1 = Label(self.frame1, text='Chromosome:', bg=DEFAULT_THEME)
        self.chrlbl1.grid(column=1, row=0)

        # Adding a Text box Entry widget for Chromosome drop-down
        chr_box = tk.StringVar()
        chr_chosen = ttk.Combobox(self.frame1, width=12, textvariable=chr_box)
        chr_chosen['values'] = ("-----Select-----", "chr1", "chr2", "chr3", "chr4", "chr5",
                                "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                                "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                                "chr18", "chr19", "chr20", "chr21", "chr22", "chrMT", "chrX", "chrY")
        chr_chosen.grid(column=2, row=0)
        chr_chosen.current(0)

        chr_chosen.bind("<<ComboboxSelected>>", MenuFunc.load_chromosome)

        self.frmChrGen = Frame(self.win, bg=DEFAULT_THEME)
        self.frmChrGen.pack(side=BOTTOM)

        self.genes_in_area = Text(self.frmChrGen, width=50, height=3, wrap=None, bg="snow")
        self.genes_in_area.grid(row=0, column=1, pady=(0, 0))
        self.genes_in_area.bind("<Enter>", Mouse.on_enter)
        self.genes_in_area.tag_configure("center", justify="center")
        self.genes_in_area.config(state=DISABLED)
        geneTxt = "G\nE\nN\nE\nS"
        self.GeneLbl = Label(self.frmChrGen, text=geneTxt, font="Helvetica 12 bold", bg=DEFAULT_THEME, fg="gray27")
        self.GeneLbl.grid(row=0, column=0, pady=(0, 0))
        coordText = "C\nO\nO\nR\nD\n"
        self.CoordLbl = Label(self.frmChrGen, text=coordText, font="Helvetica 12 bold", bg=DEFAULT_THEME, fg="gray27")
        self.CoordLbl.grid(row=0, column=2, pady=(20, 0))

        # Initialize the chromosome sequence window
        self.frameChr = Frame(self.win, bg=DEFAULT_THEME)
        self.frameChr.pack(side=BOTTOM)
        self.chrom_header = Label(self.frameChr, text="Chrom-view", font="Helvetica 18 bold", bg=DEFAULT_THEME, fg="gray27")
        self.chrom_header.grid(row=0, column=1, pady=(30, 0))
        self.chrom_area = Text(self.frameChr, width=107, height=3, wrap=NONE, inactiveselectbackground='white', bg="snow")
        self.chrom_area.grid(row=1, column=1, pady=(5, 0))
        # Add a tag to allow positions markers to indicate between bases
        self.chrom_area.tag_add("Marker", "1.0", "1.1")
        self.chrom_area.tag_config("Marker", font=("Georgia", "6", "bold"), foreground="white")
        self.chrom_area.config(state=DISABLED)
        # Chromosome area events to drag sequence
        self.chrom_area.bindtags(('Text', 'post-class-bindings', '.' 'all'))
        self.chrom_area.bind_class('post-class-bindings', '<ButtonPress-1>', Mouse.move_start)
        self.chrom_area.bind_class('post-class-bindings', '<B1-Motion>', Mouse.move_move)
        self.chrom_area.bind_class('post-class-bindings', '<ButtonRelease-1>', Mouse.on_release)
        self.chrom_area.bind_class('post-class-bindings', '<Enter>', Mouse.on_enter)

        # Initialize the Fusion area
        self.frame2 = Frame(self.win, bg=DEFAULT_THEME)
        self.frame2.pack(side=BOTTOM)
        scrollbar1 = Scrollbar(self.frame2, orient='horizontal')
        scrollbar1.grid(row=4, column=1, sticky=N + S + E + W, pady=(0, 0))
        self.fusion_header = Label(self.frame2, text="Fusion-view", font="Helvetica 18 bold", bg=DEFAULT_THEME, fg="gray27")
        self.fusion_header.grid(row=2, column=1, pady=(15, 0))
        self.fusion_area = Text(self.frame2, width=107, height=10, wrap=NONE, xscrollcommand=scrollbar1.set, bg="snow")
        scrollbar1.config(command=self.fusion_area.xview)
        self.fusion_area.grid(row=3, column=1, pady=(0, 0))

        # Makes it so that text cannot be highlight with cursor
        self.fusion_area.bindtags((str(self.fusion_area), str(self.win), "all"))
        # Add a tag to allow positions markers to indicate between bases
        self.fusion_area.tag_add("Marker", "1.0", "1.1")
        self.fusion_area.tag_config("Marker", font=("Georgia", "6", "bold"), foreground="white")
        self.fusion_area.config(state=DISABLED)

        # Fusion area events to drag sequences
        self.fusion_area.bind("<ButtonPress-1>", Mouse.move_start)
        self.fusion_area.bind("<B1-Motion>", Mouse.move_move)
        self.fusion_area.bind("<ButtonRelease-1>", Mouse.on_release)
        self.fusion_area.bind("<Enter>", Mouse.on_enter)

        self.frame3 = Frame(self.win, bg=DEFAULT_THEME)
        self.frame3.pack(side=BOTTOM)

        # Initialize the fusion boxes
        # Fusion boxes (left & right)
        self.Gene1Txt = Text(self.frame3, width=25, height=3, bg="snow")
        self.Gene1Txt.grid(row=1, column=0, padx=(0, 50), pady=(125, 0))
        self.Gene1Txt.bind("<Enter>", Mouse.on_enter)
        self.Gene1Txt.config(state=DISABLED)
        self.Gene2Txt = Text(self.frame3, width=25, height=3, bg="snow")
        self.Gene2Txt.grid(row=1, column=2, padx=(50, 0), pady=(125, 0))
        self.Gene2Txt.bind("<Enter>", Mouse.on_enter)
        self.Gene2Txt.config(state=DISABLED)

        self.Gene1Header = Label(self.frame3, text="Gene List A", fg="gray17", bg=DEFAULT_THEME, font="Helvetica 12 bold")
        self.Gene1Header.grid(row=1, column=0, padx=(0, 50), pady=(25, 0))
        self.Gene2Header = Label(self.frame3, text="Gene List B", fg="gray17", bg=DEFAULT_THEME, font="Helvetica 12 bold")
        self.Gene2Header.grid(row=1, column=2, padx=(50, 0), pady=(25, 0))

        # Initialize the summary view
        self.SummaryHeader = Label(self.frame3, text="Summary", font="Helvetica 12 bold", bg=DEFAULT_THEME, fg="gray27")
        self.SummaryHeader.grid(row=0, column=1, pady=(0, 0))
        self.SummaryTxt = Text(self.frame3, width=42, height=10, bg=DEFAULT_THEME, borderwidth=2, relief="groove")
        self.SummaryTxt.grid(row=1, column=1, pady=(0, 0))
        self.SummaryTxt.bind("<Enter>", Mouse.on_enter)
        self.SummaryTxt.config(state=DISABLED)

        self.widget_list = [self.win, self.frame1, self.refLabel, self.searchlbl1, self.chrlbl1, self.frmChrGen,
                            self.GeneLbl, self.CoordLbl, self.frameChr, self.chrom_header, self.frame2,
                            self.fusion_header, self.frame3, self.Gene1Header, self.Gene2Header,
                            self.SummaryHeader, self.SummaryTxt]

        self.update_gene_area()  # This looks in the chrom-view and updates any visible genes

        self.win.mainloop()

    # Displays the genes currently in view of the chromosome viewer
    def update_gene_area(self):
        try:
            if SEL_CHROMOSOME != "" and GENOME_LOADED:  # First check if a chromosome has been loaded
                self.genes_in_area.delete("1.0", END)  # Clear the area
                chrom_indx = self.chrom_area.index("@0,0")  # Get the first visible position in the viewer
                chrom_indx = chrom_indx.replace("1.", "3.")  # Adjust the lines to accurately represent the position
                chrom_indx = chrom_indx.replace("2.", "3.")
                chrom_seq = self.chrom_area.get(chrom_indx, END)  # Get the sequence
                chrom_seq = chrom_seq.strip(" \n")
                relative_pos = CHROMOSOME_SEQ.find(chrom_seq)  # Find the position viewer seq relative to the chromosome
                genes = BAMParser().get_gene_info(relative_pos, SEL_CHROMOSOME)
                self.genes_in_area.config(state=NORMAL)
                self.genes_in_area.delete("1.0", END)
                for gene in genes:
                    gene_string = gene.gene_name + ": " + str(gene.start) + " - " + str(gene.end)
                    self.genes_in_area.insert(INSERT, gene_string + "\n")
                self.genes_in_area.tag_add("center", "1.0", END)
                self.genes_in_area.config(state=DISABLED)
        except Exception as Ex:
            print("update_gene_area() - Error updating gene area: " + str(Ex))
            pass
        # Keep function alive so it searches for genes every 2 seconds
        self.win.after(2000, self.update_gene_area)

    # Exit GUI
    def _quit(self):
        self.win.quit()
        self.win.destroy()
        exit()


class MouseFunctions(FusionGUI):

    """
    Methods in this class deals with functions pertaining to
    a) changing mouse cursor upon entry and exit of widgets
    b) scrolling by clicking and dragging mouse in the fusion and chromosome viewers

    Functions are event based and require consistent
    reference to the GUI as it gets updated.
    """

    def __init__(self, GUI):
        self.GUI = GUI

    # GUI related functions
    @staticmethod
    def on_release(event):
        event.widget.config(cursor="")

    @staticmethod
    def on_enter(event):
        event.widget.config(cursor="")

    # Move Fusion/Chromosome viewers with mouse
    # Chromosomal viewer gets position and adjusts due to the large sequence available
    @staticmethod
    def move_start(event):  # Get the current position upon clicking in the viewer
        event.widget.scan_mark(event.x, event.y)
        event.widget.config(cursor="hand1")

    def move_move(self, event):  # Adjust the chromosome sequence when nearing the current portions beginning/end
        chrom_indx = self.GUI.chrom_area.index(INSERT)  # Default 1.0 if no chromosome is loaded
        line_num = int(str(chrom_indx)[0:1])  # Get the line number
        # Get the position in the viewer to update if reaches above or below a certain threshold
        viewer_pos = int(str(chrom_indx)[2:])
        chrom_seq = self.GUI.chrom_area.get(chrom_indx, END)
        chrom_seq = chrom_seq.strip(" \n")

        relative_pos = CHROMOSOME_SEQ.find(chrom_seq)  # Find the chromosome position only if initialized

        event.widget.scan_dragto(event.x, 0)  # This moves the viewer by dragging it

        if line_num == 3:  # Focus on the line with the sequence
            if viewer_pos <= 1000 or viewer_pos >= 9000:  # Check if viewer is reaching threshold of the loaded sequence
                # Given the current sequence, find where it is relative to the chromosome
                view_diff = relative_pos - viewer_pos
                if view_diff > 0:  # If upon reaching the threshold we still have more bases to view, load them up
                    if relative_pos >= VIEWER_SPAN:
                        # Deactivate so sequence can properly update in real-time
                        self.GUI.chrom_area.unbind_class('post-class-bindings', '<ButtonPress-1>')
                        self.GUI.chrom_area.unbind_class('post-class-bindings', '<B1-Motion>')
                        self.GUI.win.after(500, self.rebind)  # Delay rebinding to allow time for viewer to repopulate
                        # Adjust for position to keep same place
                        new_start_span = int(relative_pos - (VIEWER_SPAN/2)) + 50
                        new_end_span = int(relative_pos + (VIEWER_SPAN/2)) + 50

                        ResetFunctions(self.GUI).reset_chromosome()  # Reset chrom area and populate with new section
                        AppendFunctions(self.GUI).insert_chrom_seq(CHROMOSOME_SEQ, new_start_span, new_end_span, True)

    # Reactivates unbound functions
    def rebind(self):
        self.GUI.chrom_area.bind_class('post-class-bindings', '<ButtonPress-1>', self.move_start)
        self.GUI.chrom_area.bind_class('post-class-bindings', '<B1-Motion>', self.move_move)


class SearchFunctions(FusionGUI):

    """
    This class stores methods for finding and putting into view requested items.
    These requests are made by either searching by position or gene name in the
    search box, or when a BAM file is loaded the class will search for matching
    sequences and, if it finds any, highlight them and focus the viewer in the
    center of the requested position.
    """

    def __init__(self, GUI):
        self.GUI = GUI

    # Go to the position or gene position entered in the search bar
    def search_position(self, entered):
        entered = entered.replace(",", "")  # Format the position in case of commas
        if entered.strip() == "":
            return  # Blank search will do nothing

        if CHROMOSOME_SEQ != "":  # Check first to see if the chromosome has been loaded
            if not str(entered).isnumeric() and GENOME_LOADED:  # If reference is loaded we can find position by gene name
                entered = BAMParser().get_start_position_by_name(entered)
                entered = int(entered) + 50  # Adjust for half length of viewer
            if str(entered).isnumeric():  # If we get to here with a numeric value, go to the position
                chrom_region_start = int(entered) - 5000  # e.g. 75,758,257 - 5000 = 75,753,257
                chrom_region_end = int(entered) + 5000
                self.GUI.chrom_area.config(state=NORMAL)
                self.GUI.chrom_area.delete("1.0", END)  # Reset the viewer with the new positioned reads
                try:  # Attempt to go to entered position. If out of bounds it will alert the user.
                    AppendFunctions(self.GUI).insert_chrom_seq(CHROMOSOME_SEQ, chrom_region_start, chrom_region_end, True)
                except Exception as Ex:
                    Messages.display_msg(str(Ex))
                    pass
            elif not GENOME_LOADED:
                Messages.display_msg("Please load the reference file before searching by gene name.")

        if BAM_LOADED:  # Check if bam file has been loaded, if so, go to the position in the fusion viewer if exists
            fuse_pos, fuse_len = BAMParser().get_fusions()
            if len(fuse_pos) != 0:  # If any fusions loaded, go to the position of the topmost fusion
                adjusted_search_pos = int(entered) - int(fuse_pos[0])
                self.GUI.fusion_area.see("1." + str(adjusted_search_pos+8))

    # Looks for text and highlights it
    @staticmethod
    def search(text_widget, keyword, tag):
        pos = '1.0'
        while True:
            idx = text_widget.search(keyword, pos, END)
            if not idx:
                break
            pos = '{}+{}c'.format(idx, len(keyword))
            text_widget.tag_add(tag, idx, pos)

    # This focuses the scroll/text area on the point of fusion genes meeting
    def focus_fusion_area(self, len_diff):
        fusion_indx = "1." + str(len_diff+8)
        self.GUI.fusion_area.see(fusion_indx)


class MenuFunctions(FusionGUI):

    """
    Functions pertaining to menu items are found here.
    These are the main functions as the deal with loading
    the BAM and reference files and mapping them out.
    Additionally there are more miscellaneous functions
    such as changing the theme of the viewer, opening the
    manual and reporting errors users may come across.
    """

    def __init__(self, GUI):
        self.GUI = GUI

    # Display the document that shows how to use the application and what everything means
    @staticmethod
    def open_how_to():
        WinUse = Tk()
        WinUse.title("How to Use Application")

        frameH = Frame(WinUse)
        frameH.pack()

        helpBox = Text(frameH, width=70, height=29, wrap=WORD)
        helpBox.grid(row=0)

        HelpText1 = '''BAM ONLY: To use this application, simply select a BAM file from the 
                      'File' dropdown and the program will read it. If a fusion was found, 
                      the reads will be displayed in the Fusion Area in the middle of the 
                      application. The fusion gene will be displayed and color coded on the 
                      left side, while it's supplement will show up on the right. Positions 
                      mark the start of the base that follows and each 'marker' is 5 bases 
                      in length.'''

        HelpText = "".join(HelpText1.splitlines()) + "\n\n"

        HelpText2 = '''BAM W/ ANNOTATIONS: In order for the program to retrieve the gene 
                      names and other relevant information, a GTF reference file will need 
                      to be loaded into the program. This can be done by selecting the 
                      reference the 'Ref' tab, the type of genome, and the file associated 
                      with that genome. The program comes preloaded with Ensembles GChr38 
                      annotation file for use. You can load the reference file before 
                      loading the BAM file and vice versa. For fusions found, the Fusion 
                      gene list and the Supplement gene list will populate with discovered 
                      fusions. More information about each gene can be found by clicking 
                      'More Information' under each gene name. The summary section will 
                      show key information about the fusion that was found. 
                      ***NOTE: The reference file will have to be reloaded if loading bam 
                      files one after another to populate the gene lists and summary area.'''

        HelpText = HelpText + "".join(HelpText2.splitlines())

        HelpText3 = '''CHROMOSOMES: Each chromosome (hg38) can be loaded by using the
                      dropdown located near the top of the application. Selecting a sequence
                      will load it into the Chrom-view. Once there, the user may drag left 
                      or right and genes will populate based on what is currently in view.
                      Users may also use the search function at the top to focus on a 
                      position or go to where a specific gene is located at.'''

        HelpText = HelpText + "\n\n" + str("".join(HelpText3.splitlines()))
        helpBox.insert(INSERT, HelpText)

        WinUse.resizable(False, False)
        WinUse.mainloop()

    # TODO: Send an email to the administrator to report any bugs/areas of improvement
    @staticmethod
    def report_error():
        print("Nothing to see here.")

    # Changes the color scheme of the application
    def change_theme(self, theme):
        for wid in self.GUI.widget_list:
            wid.configure(bg=theme)

    # Given the zipped chromosome file, load the chromosome fasta into the chromosomal window
    def load_chromosome(self, event):
        global CHROM_START  # Starting position of the chromosome region. Default = 100000
        global CHROM_END  # Ending position of the chromosome region. Default = 110000
        global VIEWER_SPAN  # Number of bases that span the chromosomal viewer.
        global CHROMOSOME_SEQ
        global SEL_CHROMOSOME  # Selected chromosome

        # Distances are 10,000 apart unless search_position numbers are other than 5000 (double for viewer span)
        CHROM_START = 10000
        CHROM_END = 20000
        VIEWER_SPAN = CHROM_END - CHROM_START

        SEL_CHROMOSOME = event.widget.get().replace('chr', '')
        file = "Chromosomes/Homo_sapiens.GRCh38.dna.chromosome." + str(SEL_CHROMOSOME) + ".fa.gz"
        if "Select" not in SEL_CHROMOSOME:
            try:
                ResetFunctions(self.GUI).reset_chromosome()
                chrom_type = "CHR" + str(SEL_CHROMOSOME)  # Set the chromosome name for extraction
                init = InitFiles(file, chrom_type)
                init.create_loading()
                CHROMOSOME_SEQ = BAMParser().get_chromosome_seq()
                AppendFunctions(self.GUI).insert_chrom_seq(CHROMOSOME_SEQ, CHROM_START, CHROM_END, True)
            except Exception as Ex:
                Messages.display_error("Error loading chromosome. " + str(Ex) +
                                         "Please make sure you have the correct chromosome file: " + file)

    # Open and Load BAM File
    def load_bam_file(self):
        global BAM_LOADED
        file = filedialog.askopenfilename(filetypes=(("bam files", "*.bam"), ("All files", "*.*")))
        if file:
            try:
                ResetFunctions(self.GUI).reset_reads()  # Blank out the fusion area
                init = InitFiles(file, "SAM")
                init.create_loading()  # Create a loading bar and load the bam/sam file
                AppendFunctions(self.GUI).insert_reads(init.reads)  # Load in read information into the display
                BAM_LOADED = True
            except Exception as Ex:
                Messages.display_error("Error loading BAM file. " + str(Ex) +
                                         ". Please make sure bai file exists in the same location as the BAM.")

    # Open and Load Genome
    def load_genome_file(self, refGen):
        global FUSION_GENES
        global GENOME_LOADED

        FUSION_GENES = []

        # Ask the user to manually select a reference file
        # file = filedialog.askopenfilename(filetypes=(("gtf files", "*.gtf"), ("All files", "*.*")))
        file = "Genomes/Homo_sapiens.GRCh38.100.gtf"
        if file:
            try:
                init = InitFiles(file, "GTF")
                init.create_loading()  # Create the loading bar and load genome
                refLbl = Label(self.GUI.frameRef, text='Reference Genome ' + str(refGen) + " in Use", bg=DEFAULT_THEME)
                refLbl.grid(row=0, column=0)
                self.GUI.widget_list.append(refLbl)  # Thematic changes list
                GENOME_LOADED = True
                if BAM_LOADED:  # Check to see if a sam file has been loaded, if so reset fusion view
                    ResetFunctions(self.GUI).reset_reads()
                    BAMParser().map_fusion_reads()
            except Exception as Ex:
                Messages.display_error("Error loading GTF file. " + str(Ex) +
                                "The file may be downloaded at this address: "
                                "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz")


class Messages:

    """
    Functions to display dialog message boxes
    """

    # Notify of error
    @staticmethod
    def display_error(errMsg=""):
        messagebox.showerror("Error", errMsg)

    # Display regular message
    @staticmethod
    def display_msg(msg):
        messagebox.showinfo("Notice", msg)


class ResetFunctions(FusionGUI):

    """
    This is used to reset any viewers that are populated with read data
    including the gene information boxes.
    """

    def __init__(self, GUI):
        self.GUI = GUI

    # Clears the fusion gene area when a new file is loaded
    def reset_reads(self):
        self.GUI.fusion_area.config(state=NORMAL)
        self.GUI.Gene1Txt.config(state=NORMAL)
        self.GUI.Gene2Txt.config(state=NORMAL)
        self.GUI.SummaryTxt.config(state=NORMAL)

        self.GUI.Gene1Txt.delete('1.0', END)
        self.GUI.Gene2Txt.delete('1.0', END)
        self.GUI.SummaryTxt.delete('1.0', END)

        fuse_pos, fuse_lens = BAMParser.get_fusions()
        if len(fuse_pos) != 0:  # If fusions found we want to reset the "More information" labels here
            Gene1Label.grid_forget()
            Gene2Label.grid_forget()

        self.GUI.fusion_area.delete('1.0', END)

    # Reset the chromosome area
    def reset_chromosome(self):
        self.GUI.chrom_area.config(state=NORMAL)
        self.GUI.chrom_area.delete('1.0', END)
        self.GUI.chrom_area.config(state=DISABLED)


class AppendFunctions(FusionGUI):

    """
    The viewers are populated using functions declared in this class.
    Commonly used to shift the chromosome sequence as the viewer is
    scrolled from left to right (too many base pairs will slow the system)
    and to populate the fusion reads that were found in a BAM file.

    Class also establishes position markers for fusion and chromosome viewers
    and summary information for any found fusions as well as detailed descriptions
    of the genes themselves.
    """

    def __init__(self, GUI):
        self.GUI = GUI

    # Insert the chromosome sequence into the viewer at the bottom given a seq, start and end position
    def insert_chrom_seq(self, seq, start, end, isMiddle):
        middle_indx = "3." + str(int(VIEWER_SPAN / 2))
        if not isMiddle:
            middle_indx = "3.0"

        # The chromosome view has fixed markers every 50 bases
        # Such we must adjust the sequence to display accurately if not divisible by 50
        start = round(start/50) * 50
        remainder = (start % 50) - 1
        seq = seq[start+remainder:end+remainder]
        seq_len = len(seq)
        if seq_len == 0:
            raise Exception("Cannot get chromosome sequence: Search position out of bounds.")

        self.GUI.chrom_area.config(state=NORMAL)
        positions = BAMParser().get_position_headers(seq_len, start, 50)
        self.create_headers(seq, positions, self.GUI.chrom_area)

        # Initialize the chromosome region sequence
        self.GUI.chrom_area.insert(INSERT, "\t" + seq)
        self.GUI.chrom_area.see(middle_indx)  # Default to the middle to allow space for left and right view
        self.GUI.chrom_area.config(state=DISABLED)

    # Fusion reads are inserted line by line into the fusion area
    # Other views are created depending on loading of GTF file
    def insert_reads(self, reads):
        for read in reads:
            self.GUI.fusion_area.config(state=NORMAL)
            # Create the headers if at least one fusion was found
            if read['isFirst'] and read['fusion'] != "No Fusions Detected":
                self.create_headers(read['fusion'], read['headers'], self.GUI.fusion_area)

            if read['fusion'] == "No Fusions Detected":
                fusion = "\n\n" + read['fusion']
                self.GUI.fusion_area.tag_config('NA', foreground="indian red",
                                                font=("Georgia", "18", "bold"), justify='center')
                self.GUI.fusion_area.insert(INSERT, "\t" + fusion, 'NA')

            if read['supplement'] != "No Fusions Detected":
                self.GUI.fusion_area.insert(INSERT, "\t" + read['fusion'])
                self.GUI.fusion_area.insert(INSERT, '\n')

                len_diff = read['fusionLength'] - read['supplementLength']
                self.GUI.fusion_area.tag_config(read['fusion'], background='pale green')
                self.GUI.fusion_area.tag_config(read['supplement'], background='light salmon')
                # Highlights the fusion and supplement gene
                SearchFunctions(self.GUI).search(self.GUI.fusion_area, read['fusion'], read['fusion'])
                SearchFunctions(self.GUI).search(self.GUI.fusion_area, read['supplement'], read['supplement'])

                SearchFunctions(self.GUI).focus_fusion_area(len_diff)

            # Lock the text box
            self.GUI.fusion_area.config(state=DISABLED)
            if read['fusionName'] != "" and read['fusionName'] not in FUSION_GENES:
                self.add_gene_info(read['fusionName'], read['supplementName'], read['fusionId'], read['supplement'])
                self.add_fusion_summary(read['fusionName'], read['supplementName'], read['fusionPosition'],
                                        read['supplementPosition'], read['fusionLength'], read['supplementLength'])
                FUSION_GENES.append(read['fusionName'])
                FUSION_GENES.append(read['supplementName'])

    # Place the position markers for the first row read
    def create_headers(self, seq, positions, viewer):
        length = len(seq)
        dist_between = 9  # Number of markers before a positional marker is set
        pos_marker = "\u2502"
        mini_markers = "    \u2758" * dist_between
        viewer.tag_config('Header', background="gray81")
        viewer.insert(INSERT, "|", "Marker")  # Small pipe shifts position markers slightly so to land in between bases
        viewer.insert(INSERT, " " * length, 'Header')  # Initialize the header region so indeces may be referenced
        viewer.insert(INSERT, "\n")
        viewer.insert(INSERT, "|", "Marker")  # Twice for the positions and the line markers themselves
        viewer.insert(INSERT, " " * length, 'Header')
        viewer.insert(INSERT, "\n")

        start_position = positions[0]

        for x in range(len(positions)):
            position_chars = len(str(positions[x]))  # Get the number of characters in the position (as string)
            pos_adjusted_indx = round(position_chars/2)
            adjusted_start_indx = 8 - pos_adjusted_indx  # Shifts position into the center of marker (8 = length of tab)

            if x == 0:
                viewer.insert("1." + str(adjusted_start_indx), positions[x])
                viewer.insert("2.8", pos_marker + mini_markers)
            elif x != len(positions):
                position_indx = "1." + str((positions[x] - start_position) + adjusted_start_indx)  # Insert the position
                viewer.insert(position_indx, positions[x])
                position_indx = "2." + str((positions[x] - start_position) + 8)  # Insert the marker
                if x != len(positions)-1:
                    viewer.insert(position_indx, pos_marker + mini_markers)  # Add mini-markers in between positions
                else:
                    viewer.insert(position_indx, pos_marker)
        viewer.delete("1." + str(length+20), '1.99999999999')  # Remove extra spaces and leave room for the end position
        viewer.delete("2." + str(length+20), '2.99999999999')

    # List the fusion genes and add additional info
    def add_gene_info(self, Gene1Name, Gene2Name, Gene1ID, Gene2ID):
        global Gene1Label
        global Gene2Label

        Gene1Link = "https://uswest.ensembl.org/Homo_sapiens/Gene/Summary?g=" + str(Gene1ID)
        Gene2Link = "https://uswest.ensembl.org/Homo_sapiens/Gene/Summary?g=" + str(Gene2ID)

        self.GUI.Gene1Txt.config(state=NORMAL)
        self.GUI.Gene1Txt.insert(INSERT, Gene1Name)
        self.GUI.Gene1Txt.tag_add("here", "1.0", "end")
        self.GUI.Gene1Txt.tag_config("here", background="pale green", foreground="black")
        self.GUI.Gene1Txt.tag_configure("center", justify="center")
        self.GUI.Gene1Txt.tag_add("center", "1.0", "end")
        self.GUI.Gene1Txt.insert(INSERT, "\n\n")
        self.GUI.Gene1Txt.config(state=DISABLED)

        var1 = StringVar()
        Gene1Label = Label(self.GUI.frame3, width=24, height=1, textvariable=var1,
                           bg="pale green", fg="RoyalBlue3", font="Helvetica 11 bold")
        Gene1Label.grid(row=1, column=0, padx=(0, 50), pady=(150, 0))
        var1.set("More information")
        Gene1Label.bind("<Button-1>", lambda e: AdditionalFunctions(self.GUI).show_more_info(Gene1Name, Gene1Link))

        self.GUI.Gene2Txt.config(state=NORMAL)
        self.GUI.Gene2Txt.insert(INSERT, Gene2Name)
        self.GUI.Gene2Txt.tag_add("here", "1.0", "end")
        self.GUI.Gene2Txt.tag_config("here", background="light salmon", foreground="black")
        self.GUI.Gene2Txt.tag_configure("center", justify="center")
        self.GUI.Gene2Txt.tag_add("center", "1.0", "end")
        self.GUI.Gene2Txt.insert(INSERT, "\n\n")
        self.GUI.Gene2Txt.config(state=DISABLED)

        var2 = StringVar()
        Gene2Label = Label(self.GUI.frame3, width=24, height=1, textvariable=var2,
                           bg="light salmon", fg="RoyalBlue3", font="Helvetica 11 bold")
        Gene2Label.grid(row=1, column=2, padx=(50, 0), pady=(150, 0))
        var2.set("More information")
        Gene2Label.bind("<Button-1>", lambda e: AdditionalFunctions(self.GUI).show_more_info(Gene2Name, Gene2Link))

    # If fusion is found, display the summary
    def add_fusion_summary(self, Gene1Name, Gene2Name, Gene1Pos, Gene2Pos, Gene1Length, Gene2Length):
        len_diff = Gene1Length - Gene2Length
        intersection = Gene1Pos + len_diff
        Gene1End = str(Gene1Pos + Gene1Length)
        Gene1Chrom = str(BAMParser().get_chromosome_by_name(Gene1Name[0]))
        Gene2End = str(Gene2Pos + Gene2Length)
        Gene2Chrom = str(BAMParser().get_chromosome_by_name(Gene2Name[0]))

        self.GUI.SummaryTxt.config(state=NORMAL)
        self.GUI.SummaryTxt.insert(INSERT, Gene1Name[0] + "-" + Gene2Name[0] + "\n\n")
        self.GUI.SummaryTxt.insert(INSERT, "Gene A Read Coordinates:\n")
        self.GUI.SummaryTxt.insert(INSERT, "Chromosome " + str(Gene1Chrom) + ": " +
                                   str(Gene1Pos) + " - " + str(Gene1End) + "\n\n")
        self.GUI.SummaryTxt.insert(INSERT, "Fusion Intersection: " + str(intersection) + "\n\n")
        self.GUI.SummaryTxt.insert(INSERT, "Gene B Read Coordinates:\n")
        self.GUI.SummaryTxt.insert(INSERT, "Chromosome " + str(Gene2Chrom) + ": " +
                                   str(Gene2Pos) + " - " + str(Gene2End) + "\n\n")
        self.GUI.SummaryTxt.tag_configure("center", justify="center")
        self.GUI.SummaryTxt.tag_add("center", "1.0", "end")
        self.GUI.SummaryTxt.config(state=DISABLED)


class AdditionalFunctions(FusionGUI):
    """
    Miscellaneous functions dealing with GUI features not relevant
    to any of the above classes.
    """
    def __init__(self, GUI):
        self.GUI = GUI

    # Open a window with extra information about the specific gene
    @staticmethod
    def show_more_info(GeneName, GeneLink):
        subWin = Tk()
        subWin.title(GeneName)

        frame = Frame(subWin)
        frame.pack()

        LinkLabel = Label(frame, width=45, height=1, text="Click for Online Resources", background="light cyan",
                          foreground="blue", font="Helvetica 11 bold")
        LinkLabel.bind("<Button-1>", lambda e: webbrowser.open_new(GeneLink))
        LinkLabel.pack(side=TOP)

        infoBox = Text(frame,  width=45, height=15)

        GeneStart, GeneEnd = BAMParser().get_positions_by_name(GeneName[0])

        infoBox.insert(INSERT, "ID: " + str(BAMParser().get_id_by_name(GeneName[0])) + "\n")
        infoBox.insert(INSERT, "Positions: " + str(GeneStart) + " - " + str(GeneEnd) + "\n")
        infoBox.insert(INSERT, "Chromosome: " + str(BAMParser().get_chromosome_by_name(GeneName[0])) + "\n")

        # Exon regions are identified bases on transcripts
        infoBox.insert(INSERT, "Transcripts: " + "\n")
        trans_ids = BAMParser().get_transcript_ids_by_name(GeneName[0])
        for t_id in trans_ids:
            t_start, t_end = BAMParser().get_transcript_positions_by_id(t_id)
            t_name = BAMParser().get_transcript_name_by_id(t_id)
            infoBox.insert(INSERT, "  " + str(t_name) + "(" + t_id + ")\n  Positions: "
                           + str(t_start) + " - " + str(t_end) + "\n")
            infoBox.insert(INSERT, "    Exon/Intron Boundaries:\n")
            exon_start, exon_end = BAMParser().get_exon_regions_by_transcript_id(t_id)
            for i in range(len(exon_start)):
                # If the start position does not equal start position of exon, insert intron
                if int(t_start) != exon_start[0] and i == 0:
                    infoBox.insert(INSERT, "      I: " + str(t_start) + " - " + str(exon_start[i]-1) + "\n")
                infoBox.insert(INSERT, "      E: " + str(exon_start[i]) + " - " + str(exon_end[i]) + "\n")
                if i != len(exon_start)-1:  # If not on the last exon, look ahead to get intron position
                    # Insert intron positions as before the start of the next exon and after the end of the current one
                    infoBox.insert(INSERT, "      I: " + str(exon_end[i]+1) + " - " + str(exon_start[i+1]-1) + "\n")

        infoBox.pack(side=BOTTOM)
        infoBox.config(state=DISABLED)

        subWin.update()
        subWin.attributes("-topmost", True)
        subWin.resizable(False, False)
        subWin.mainloop()


if __name__ == "__main__":
    FusionGUI().runblue()
