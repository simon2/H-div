block:
	a=1 ; while [[ $$a -le 20 ]] ; do qsub jobscript.sh ; ((a = a + 1)) ; done

clean:
	rm -f B* core.* *~

omp:
	b=1 ; while [[ $$b -le 20 ]] ; do qsub omp.sh ; ((b = b + 1)) ; done

resb:
	cat B*.o* | grep "block"

resc:
	cat B*.o* | grep "cluster" | grep -v "block"

seq:
	c=1 ; while [[ $$c -le 20 ]] ; do qsub seq.sh ; ((c = c + 1)) ; done

comb:
	icc -std=c99 -O3 -o hdcm hacapk_division_cilk_merge.c -lm

como:
	icc -std=c99 -O3 -qopenmp -o hdo hacapk_division_omp.c -lm

coms:
	icc -std=c99 -O3 -o hd hacapk_division.c -lm

tssrun:
	tssrun -A p=1:c=36:t=36:m=120G ./hdcm 36

ress:
	cat B*.o* | grep "sequtial time spent"

resj:
	cat B*.o* | grep "jump access"

resu:
	cat B*.o* | grep "succession access"
