# Makefile

# Define targets and dependencies
TEX_FILE = paper.tex
POSTER_TEX = poster.tex
BIB_FILE = bibliography.bib
DOCX_FILE = paper.docx
PDF_FILE = paper.pdf
POSTER_FILE = poster.pdf
IMAGES_DIR = images/
ALL_IMAGES = $(wildcard $(IMAGES_DIR)/*)

all: $(PDF_FILE) $(DOCX_FILE)

paper: $(PDF_FILE) $(DOCX_FILE)

poster: $(POSTER_FILE)

$(POSTER_FILE): $(POSTER_TEX) $(BIB_FILE) $(ALL_IMAGES)
	latexmk $(POSTER_TEX)

# Target: paper.pdf
$(PDF_FILE): $(TEX_FILE) $(BIB_FILE) $(ALL_IMAGES)
	latexmk -pdf $(TEX_FILE)

# Target: paper.docx
$(DOCX_FILE): $(TEX_FILE) $(BIB_FILE) $(IMAGES_DIR)
	pandoc $(TEX_FILE) -o $(DOCX_FILE)

# Phony target to clean up temporary files
clean:
	latexmk -C
	rm *.bbl
	rm *.run.xml
	rm *.fdb_latexmk
	rm *.aux
	rm *.fls
	rm *.log
	rm -f $(DOCX_FILE)
