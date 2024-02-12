# Makefile

# Define targets and dependencies
TEX_FILE = paper.tex
BIB_FILE = bibliography.bib
DOCX_FILE = paper.docx
PDF_FILE = paper.pdf
IMAGES_DIR = images/

all: $(PDF_FILE) $(DOCX_FILE)

# Target: paper.pdf
$(PDF_FILE): $(TEX_FILE) $(BIB_FILE) $(IMAGES_DIR)
	latexmk -pdf $(TEX_FILE)

# Target: paper.docx
$(DOCX_FILE): $(TEX_FILE) $(BIB_FILE) $(IMAGES_DIR)
	pandoc $(TEX_FILE) -o $(DOCX_FILE)

# Phony target to clean up temporary files
clean:
	latexmk -c
	rm -f $(DOCX_FILE)
