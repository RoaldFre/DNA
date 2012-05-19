TARGET=report

all: tex

ALL: 
	make -B all -i
	make -B bib 
	make -B all -i
	make -B all 

tex:
	pdflatex ${TARGET}.tex
	rm ${TARGET}.log

bib:
	bibtex ${TARGET}

clean:
	rm -f ${TARGET}.{log,aux,bbl,blg}

