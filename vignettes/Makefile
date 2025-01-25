all:	vign

vign:	tr_2004_01.tex
	pdflatex tr_2004_01
	pdflatex tr_2004_01
	make clean

pdf:	tr_2004_01.tex
	echo "Sweave('tr_2004_01.Rnw')" | nohup R-rel --vanilla > out.txt
	pdflatex tr_2004_01
	pdflatex tr_2004_01
	make clean
	acroread tr_2004_01.pdf &

clean:
	rm -f *.out *.bbl *.log *.aux *.blg *.toc  *~
	rm -f *.dvi *.ps
	rm -f tr_2004_01-* out.txt Rplots.ps Rplots.pdf
