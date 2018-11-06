(c-exp "#include <stdio.h>")
(c-exp "#include <stdlib.h>")
(c-exp "#include <math.h>")
(c-exp "#include <time.h>")
(c-exp "#include <sys/time.h>")
(%defconstant INPUT_DEFAULT "bem_data/input_50ms.txt")

;; *********define cluster************
(deftype cluster (struct cluster))
(def (struct cluster)
  (decl ndim int)
  (decl nstrt int)
  (decl nsize int)
  (decl ndpth int)
  (decl nnson int)
  (decl nmbr int) 
  (decl ndscd int)
  (decl bmin (ptr double))
  (decl bmax (ptr double))
  (decl zwdth double)
  (decl pc-sons (ptr (ptr cluster))))

(deftype leafmtx (struct leafmtx))
(def (struct leafmtx) 
  (decl ltmtx int)
  (decl kt int)
  (decl nstrtl int) 
  (decl ndl int)
  (decl nstrtt int) 
  (decl ndt int) 
  (decl a1 (ptr double))
  (decl a2 (ptr double)))

(deftype leafmtxp (struct leafmtxp))
(def (struct leafmtxp)
    (decl nlf int)
  (decl nlfkt int))

(decl (supermatrix-construction-cog-leafmtrx st-leafmtxp gmid param lod lnmtx nofc nffc ndim)
    (fn void (ptr leafmtxp) (ptr (array double 3)) (array double) (ptr int) (ptr int) int int int))
(decl (qsort-col-leafmtx st-leafmtx first last)
    (fn void (ptr leafmtx) int int))
(decl (qsort-row-leafmtx st-leafmtx first last)
    (fn void (ptr leafmtx) int int))
(decl (med3 nl nr nlr2)
    (fn int int int int))
(decl (create-leafmtx st-leafmtx stcltl st-cltt param lnmtx nffc nlf)
    (fn void (ptr leafmtx) (ptr cluster) (ptr cluster) (array double) (ptr int) int (ptr int)))
(decl (dist-2cluster st-cltl stcltt)
    (fn double (ptr cluster) (ptr cluster)))
(decl (count-lntmx st-cltl st-cltt param lnmtx nffc)
    (fn void (ptr cluster) (ptr cluster) (array double) (ptr int) int))
(decl (cal-bndbox-cog st-clt zgmid lod nofc)
    (fn void (ptr cluster) (ptr (array double 3)) (ptr int) int))
(decl (set-bndbox-cog st-clt zgmid lod nofc)
    (fn void (ptr cluster) (ptr (array double 3)) (ptr int) int))
(decl (create-cluster nmbr ndpth nstrt nsize ndim nson)
    (fn (ptr cluster) int int int int int int))
(decl (free-st-clt st-clt)
    (fn void (ptr cluster)))
(decl (create-ctree-ssgeom st-clt zgmid param lod ndpth ndscd nsrt nd md ndim nclst)
    (fn (ptr cluster) (ptr cluster) (ptr (array double 3)) (array double) (ptr int) int int int int int int int))
(decl (csym::get-wall-time) (fn double))
(decl (csym::get-cpu-time) (fn double))
(decl (checkClusterTree f st-clt) (fn void (ptr FILE) (ptr cluster)))

(decl depth-max int)
(decl count-node int)

(def (main argc argv) (fn int int (ptr (ptr char)))
  (decl fname (ptr char))
  (decl file (ptr FILE))
  (def countOfNode int 0)
  (def count int 0)
  (decl i int)
  (decl coordOfNode (ptr (array double 3)))
  (decl coordOfFace (ptr (array double 3)))
  (= fname (if-exp (>= argc 2) (aref argv 1) INPUT-DEFAULT))
  (= file (fopen fname "r"))
  (if (== file NULL)
      (begin
        (csym::printf "Error: Unable to input file '%s'!~%" fname)
        (csym::exit 99))
      (begin
        (decl line (array char 100))
        (csym::fgets line (sizeof line) file)
        (csym::sscanf line "%d" (ptr countOfNode))
        (= coordOfNode (cast (ptr (array double 3))
                         (csym::malloc (* (* countOfNode 3) (sizeof double)))))
        (for ((= i 0) (< i countOfNode) (inc i))
          (csym::fgets line (sizeof line) file)
          (decl x double)
          (decl y double)
          (decl z double)
          (csym::sscanf line "%lf %lf %lf" (ptr x) (ptr y) (ptr z))
          (= (aref (aref coordOfNode i) 0) x)
          (= (aref (aref coordOfNode i) 1) y)
          (= (aref (aref coordOfNode i) 2) z))
        (csym::fgets line (sizeof line) file)
        (csym::sscanf line "%d" (ptr count))
        (= coordOfFace (cast (ptr (array double 3))
                         (csym::malloc (* (* count 3) (sizeof double)))))
        (csym::fgets line (sizeof line) file)
        (csym::fgets line (sizeof line) file)
        (csym::fgets line (sizeof line) file)
        (csym::printf "count:%d~%" count)
        (for ((= i 0) (< i count) (inc i))
          (csym::fgets line (sizeof line) file)
          (decl X int) (decl Y int) (decl Z int)
          (csym::sscanf line "%d %d %d" (ptr X) (ptr Y) (ptr Z))
          (= (aref (aref coordOfFace i) 0)
             (/ (+ (+ (aref (aref coordOfNode X) 0) (aref (aref coordOfNode Y) 0))
                   (aref (aref coordOfNode Z) 0))
                3))
          (= (aref (aref coordOfFace i) 1)
             (/ (+ (+ (aref (aref coordOfNode X) 1) (aref (aref coordOfNode Y) 1))
                   (aref (aref coordOfNode Z) 1))
                3))
          (= (aref (aref coordOfFace i) 2)
             (/ (+ (+ (aref (aref coordOfNode X) 2) (aref (aref coordOfNode Y) 2))
                   (aref (aref coordOfNode Z) 2))
                3)))))
  (csym::fclose file)
  (csym::free coordOfNode)
  (decl param (array double 100))
  (for ((= i 0) (< i 100) (inc i))
    (= (aref param i) 0))
  (= (aref param 21) 10.0)
  (= (aref param 31) 1.1)
  (= (aref param 41) 15.0)
  (= (aref param 51) 2.0)
  (decl st-leafmtxp (ptr leafmtxp))
  (decl lod (ptr int))
  (decl lnmtx (ptr int))
  (def nofc int count)
  (def nffc int 1)
  (def ndim int 3)
  (= lnmtx (cast (ptr int) (csym::malloc (* 3 (sizeof int)))))
  (for ((= i 0) (< i 3) (inc i))
    (= (aref lnmtx i) 0))
  (= st-leafmtxp (cast (ptr leafmtxp) (csym::malloc (sizeof leafmtxp))))
  (= lod (cast (ptr int) (csym::malloc (* nofc (sizeof int)))))
  (for ((= i 0) (< i nofc) (inc i))
    (= (aref lod i) 0))
  (supermatrix-construction-cog-leafmtrx st-leafmtxp coordOfFace param lod lnmtx
                                         nofc nffc ndim)
  (return 0))

(def (supermatrix-construction-cog-leafmtrx st-leafmtxp gmid param lod lnmtx nofc nffc ndim)
    (fn void (ptr leafmtxp) (ptr (array double 3)) (array double) (ptr int) (ptr int) int int int)
  (def st-clt (ptr cluster) (cast (ptr cluster) (csym::malloc (sizeof cluster))))
  (decl i int)
  (decl nfl int)
  (decl nflkt int)
  (decl ip int) 
  (decl il int)
  (decl ig int)
  (def nd int (* nofc nffc))
  (decl lodfc (ptr int))
  (decl st-leafmtx (ptr leafmtx))
  (= lodfc (cast (ptr int) (csym::malloc (* nofc (sizeof int)))))
  (for ((= il 0) (< il nofc) (inc il))
    (= (aref lodfc il) il))
  (def nsrt int 0)
  (def ndf int nofc)
  (def nclst int 0)
  (def ndpth int 0)
  (def ndscd int 0)
  (= depth-max 0)
  (= count-node 0)
  (decl start double)
  (decl end double)
  (decl spent double)
  (= start (csym::get-wall-time))
  (= st-clt
     (create-ctree-ssgeom st-clt gmid param lodfc ndpth ndscd nsrt
                          ndf nofc ndim nclst))
  (= end (csym::get-wall-time))
  (= spent (- end start))
  (csym::printf "cluster tree time spent:%.10f~%" spent)
  (set-bndbox-cog st-clt gmid lodfc nofc)
  (= ndpth 0)
  (= start (csym::get-wall-time))
  (count-lntmx st-clt st-clt param lnmtx nffc)
  (= end (csym::get-wall-time))
  (= spent (- end start))
  (csym::printf "count time:%.10f~%" spent)
  (= (fref (mref st-leafmtxp) nlfkt) (aref lnmtx 0))
  (def nlf int (+ (aref lnmtx 0) (aref lnmtx 1)))
  (= st-leafmtx (cast (ptr leafmtx) (csym::malloc (* nlf (sizeof leafmtx)))))
  (= (fref (mref st-leafmtxp) nlf) nlf)
  (csym::printf "nlf:%d~%" nlf)
  (= nlf 0)
  (= start (csym::get-wall-time))
  (create-leafmtx st-leafmtx st-clt st-clt param lnmtx nffc (ptr nlf))
  (= end (csym::get-wall-time))
  (= spent (- end start))
  (csym::printf "nlf:%d~%" nlf)
  (csym::printf "block cluster tree time spent:%.10f~%" spent)
  (csym::printf "depth_max:%d  count_node:%d~%" depth-max count-node))

(def (qsort-row-leafmtx st-leafmtx first last)
    (fn void (ptr leafmtx) int int)
  (decl i int) 
  (decl j int)
  (decl pivot int)
  (decl st-www leafmtx)
  (if (< first last)
      (begin 
        (= pivot first)
        (= i first)
        (= j last)
        (while (< i j)
          (while (and (<= (fref (aref st-leafmtx i) nstrtl)
                          (fref (aref st-leafmtx pivot) nstrtl))
                      (< i last))
            (inc i))
          (while (> (fref (aref st-leafmtx j) nstrtl)
                    (fref (aref st-leafmtx pivot) nstrtl))
            (dec j))
          (if (< i j)
              (begin
                (= st-www (aref st-leafmtx i))
                (= (aref st-leafmtx i) (aref st-leafmtx j))
                (= (aref st-leafmtx j) st-www))))
        (= st-www (aref st-leafmtx pivot))
        (= (aref st-leafmtx pivot) (aref st-leafmtx j))
        (= (aref st-leafmtx j) st-www)
        (qsort-row-leafmtx st-leafmtx first (- j 1))
        (qsort-row-leafmtx st-leafmtx (+ j 1) last))))

(def (qsort-col-leafmtx st-leafmtx first last) (fn void (ptr leafmtx) int int)
  (decl i int)
  (decl j int)
  (decl pivot int)
  (decl st-www leafmtx)
  (if (< first last)
      (begin 
        (= pivot first)
        (= i first)
        (= j last)
        (while (< i j)
          (while
              (and (<= (fref (aref st-leafmtx i) nstrtt)
                       (fref (aref st-leafmtx pivot) nstrtt))
                   (< i last))
            (inc i))
          (while
              (> (fref (aref st-leafmtx j) nstrtt) (fref (aref st-leafmtx pivot) nstrtt))
        (dec j))
       (if (< i j)
           (begin (= st-www (aref st-leafmtx i))
            (= (aref st-leafmtx i) (aref st-leafmtx j))
            (= (aref st-leafmtx j) st-www))))
      (= st-www (aref st-leafmtx pivot))
      (= (aref st-leafmtx pivot) (aref st-leafmtx j))
      (= (aref st-leafmtx j) st-www)
      (qsort-col-leafmtx st-leafmtx first (- j 1))
      (qsort-col-leafmtx st-leafmtx (+ j 1) last))))

(def (med3 nl nr nlr2) (fn int int int int)
  (decl med3 int)
  (if (< nl nr)
      (begin
        (if (< nr nlr2)
            (begin (= med3 nr))
            (if (< nlr2 nl)
                (begin (= med3 nl))
                (begin (= med3 nlr2)))))
      (begin
        (if (< nlr2 nr)
            (begin (= med3 nr))
            (if (< nl nlr2) 
                (begin (= med3 nl))
                (begin (= med3 nlr2))))))
  (return med3))

(def (create-leafmtx st-leafmtx st-cltl st-cltt param lnmtx nffc nlf)
    (fn void (ptr leafmtx) (ptr cluster) (ptr cluster) (array double) (ptr int) int (ptr int))
  (def ndl int (* (fref (mref st-cltl) nsize) nffc))
  (def ndt int (* (fref (mref st-cltt) nsize) nffc))
  (def nstrtl int (fref (mref st-cltl) nstrt))
  (def nstrtt int (fref (mref st-cltt) nstrt))
  (def nnsonl int (fref (mref st-cltl) nnson))
  (def nnsont int (fref (mref st-cltt) nnson))
  (decl il int) (decl it int)
  (def nleaf double (aref param 41))
  (def zeta double (aref param 51))
  (def zdistlt double (dist-2cluster st-cltl st-cltt))
  (if (and (or (<= (* (fref (mref st-cltl) zwdth) zeta) zdistlt)
               (<= (* (fref (mref st-cltt) zwdth) zeta) zdistlt))
           (and (>= ndl nleaf) (>= ndt nleaf)))
      (begin
        (= (fref (aref st-leafmtx (mref nlf)) nstrtl) nstrtl)
        (= (fref (aref st-leafmtx (mref nlf)) ndl) ndl)
        (= (fref (aref st-leafmtx (mref nlf)) nstrtt) nstrtt)
        (= (fref (aref st-leafmtx (mref nlf)) ndt) ndt)
        (= (fref (aref st-leafmtx (mref nlf)) kt) 0)
        (= (fref (aref st-leafmtx (mref nlf)) ltmtx) 1)
        (= (mref nlf) (+ (mref nlf) 1)))
      (begin
        (if (or (== nnsonl 0) (== nnsont 0) (<= ndl nleaf) (<= ndt nleaf))
            (begin
              (= (fref (aref st-leafmtx (mref nlf)) nstrtl) nstrtl)
              (= (fref (aref st-leafmtx (mref nlf)) ndl) ndl)
              (= (fref (aref st-leafmtx (mref nlf)) nstrtt) nstrtt)
              (= (fref (aref st-leafmtx (mref nlf)) ndt) ndt)
              (= (fref (aref st-leafmtx (mref nlf)) ltmtx) 2)
              (= (mref nlf) (+ (mref nlf) 1)))
            (begin
              (for ((= il 0) (< il nnsonl) (inc il))
                (for ((= it 0) (< it nnsont) (inc it))
                  (create-leafmtx st-leafmtx (aref (fref (mref st-cltl) pc-sons) il)
                                  (aref (fref (mref st-cltt) pc-sons) it) param lnmtx nffc
                                  nlf))))))))

(def (dist-2cluster st-cltl st-cltt) (fn double (ptr cluster) (ptr cluster))
 (def zs double 0.0) (decl id int)
 (for ((= id 0) (< id (fref (mref st-cltl) ndim)) (inc id))
   (if (< (aref (fref (mref st-cltl) bmax) id) (aref (fref (mref st-cltt) bmin) id))
       (begin
         (= zs
            (+ zs
               (* (- (aref (fref (mref st-cltt) bmin) id)
                     (aref (fref (mref st-cltl) bmax) id))
                  (- (aref (fref (mref st-cltt) bmin) id)
                     (aref (fref (mref st-cltl) bmax) id))))))
       (if (< (aref (fref (mref st-cltt) bmax) id)
              (aref (fref (mref st-cltl) bmin) id))
           (begin
             (= zs
                (+ zs
                   (* (- (aref (fref (mref st-cltl) bmin) id)
                         (aref (fref (mref st-cltt) bmax) id))
                      (- (aref (fref (mref st-cltl) bmin) id)
                         (aref (fref (mref st-cltt) bmax) id)))))))))
 (return (sqrt zs)))

(def (count-lntmx st-cltl st-cltt param lnmtx nffc)
    (fn void (ptr cluster) (ptr cluster) (array double) (ptr int) int)
  (decl il int)
  (decl it int)
  (def ndl int (* (fref (mref st-cltl) nsize) nffc))
  (def ndt int (* (fref (mref st-cltt) nsize) nffc))
  (def nstrtl int (fref (mref st-cltl) nstrt))
  (def nstrtt int (fref (mref st-cltt) nstrt))
  (def nnsonl int (fref (mref st-cltl) nnson))
  (def nnsont int (fref (mref st-cltt) nnson))
  (def nleaf double (aref param 41))
  (def zeta double (aref param 51))
  (def zdistlt double (dist-2cluster st-cltl st-cltt))
  (if (and (or (<= (* (fref (mref st-cltl) zwdth) zeta) zdistlt)
               (<= (* (fref (mref st-cltt) zwdth) zeta) zdistlt))
           (and (>= ndl nleaf) (>= ndt nleaf)))
      (begin
        (= (aref lnmtx 0) (+ (aref lnmtx 0) 1)))
      (begin
        (if (or (== nnsonl 0) (== nnsont 0) (<= ndl nleaf) (<= ndt nleaf))
            (begin
              (= (aref lnmtx 1) (+ (aref lnmtx 1) 1)))
            (begin
              (= (aref lnmtx 2) (+ (aref lnmtx 2) 1))
              (for ((= il 0) (< il nnsonl) (inc il))
                (for ((= it 0) (< it nnsont) (inc it))
                  (count-lntmx (aref (fref (mref st-cltl) pc-sons) il)
                               (aref (fref (mref st-cltt) pc-sons) it) param lnmtx nffc))))))))

(def (cal-bndbox-cog st-clt zgmid lod nofc)
    (fn void (ptr cluster) (ptr (array double 3)) (ptr int) int)
  (def ndim int (fref (mref st-clt) ndim))
  (decl id int)
  (decl il int)
  (= (fref (mref st-clt) bmin) (cast (ptr double) (csym::malloc (* 3 (sizeof double)))))
  (= (fref (mref st-clt) bmax) (cast (ptr double) (csym::malloc (* 3 (sizeof double)))))
  (def zeps double 1.0e-5)
  (if (> (fref (mref st-clt) nnson) 0)
      (begin
        (for ((= id 0) (< id ndim) (inc id))
          (= (aref (fref (mref st-clt) bmin) id)
             (aref (fref (mref (aref (fref (mref st-clt) pc-sons) 0)) bmin) id))
          (= (aref (fref (mref st-clt) bmax) id)
             (aref (fref (mref (aref (fref (mref st-clt) pc-sons) 0)) bmax) id)))
        (for ((= il 1) (< il (fref (mref st-clt) nnson)) (inc il))
          (for ((= id 0) (< id ndim) (inc id))
            (if (< (aref (fref (mref (aref (fref (mref st-clt) pc-sons) il)) bmin) id)
                   (aref (fref (mref st-clt) bmin) id))
                (begin
                  (= (aref (fref (mref st-clt) bmin) id)
                     (aref (fref (mref (aref (fref (mref st-clt) pc-sons) il)) bmax) id))))
            (if (< (aref (fref (mref st-clt) bmax) id)
                   (aref (fref (mref (aref (fref (mref st-clt) pc-sons) il)) bmax) id))
                (begin
                  (= (aref (fref (mref st-clt) bmax) id)
                     (aref (fref (mref (aref (fref (mref st-clt) pc-sons) il)) bmax) id)))))))
      (begin
        (for ((= id 0) (< id ndim) (inc id))
          (= (aref (fref (mref st-clt) bmin) id) (aref (aref zgmid (aref lod 0)) id))
          (= (aref (fref (mref st-clt) bmax) id) (aref (aref zgmid (aref lod 0)) id)))
        (for ((= id 0) (< id ndim) (inc id))
          (for ((= il 1) (< il (fref (mref st-clt) nsize)) (inc il))
            (if (< (aref (aref zgmid (aref lod il)) id)
                   (aref (fref (mref st-clt) bmin) id))
                (begin
                  (= (aref (fref (mref st-clt) bmin) id)
                     (aref (aref zgmid (aref lod il)) id))))
            (if (< (aref (fref (mref st-clt) bmax) id)
                   (aref (aref zgmid (aref lod il)) id))
                (begin
                  (= (aref (fref (mref st-clt) bmax) id)
                     (aref (aref zgmid (aref lod il)) id))))))))
  (def zwdth double
    (* (- (aref (fref (mref st-clt) bmax) 0) (aref (fref (mref st-clt) bmin) 0))
       (- (aref (fref (mref st-clt) bmax) 0) (aref (fref (mref st-clt) bmin) 0))))
  (for ((= id 1) (< id ndim) (inc id))
    (= zwdth
       (+ zwdth
          (* (- (aref (fref (mref st-clt) bmax) id)
                (aref (fref (mref st-clt) bmin) id))
             (- (aref (fref (mref st-clt) bmax) id)
                (aref (fref (mref st-clt) bmin) id))))))
  (= zwdth (sqrt zwdth))
  (for ((= id 0) (< id ndim) (inc id))
    (def bdiff double
      (- (aref (fref (mref st-clt) bmax) id) (aref (fref (mref st-clt) bmin) id)))
    (if (< bdiff (* zeps zwdth))
        (begin
          (= (aref (fref (mref st-clt) bmax) id)
             (+ (aref (fref (mref st-clt) bmax) id)
                (* 0.5 (- (* zeps zwdth) bdiff))))
          (= (aref (fref (mref st-clt) bmin) id)
             (- (aref (fref (mref st-clt) bmin) id)
                (* 0.5 (- (* zeps zwdth) bdiff)))))))
  (= zwdth
     (* (- (aref (fref (mref st-clt) bmax) 0) (aref (fref (mref st-clt) bmin) 0))
        (- (aref (fref (mref st-clt) bmax) 0) (aref (fref (mref st-clt) bmin) 0))))
  (for ((= id 1) (< id ndim) (inc id))
    (= zwdth
       (* (- (aref (fref (mref st-clt) bmax) id) (aref (fref (mref st-clt) bmin) id))
          (- (aref (fref (mref st-clt) bmax) id)
             (aref (fref (mref st-clt) bmin) id)))))
  (= (fref (mref st-clt) zwdth) (sqrt zwdth)))

(def (set-bndbox-cog st-clt zgmid lod nofc)
    (fn void (ptr cluster) (ptr (array double 3)) (ptr int) int) (decl ic int) (decl l int)
    (for ((= ic 0) (< ic (fref (mref st-clt) nnson)) (inc ic))
      (if (== ic 0) (begin (= l 0))
          (begin
            (= l
               (+ l (fref (mref (aref (fref (mref st-clt) pc-sons) (- ic 1))) nsize)))))
      (set-bndbox-cog (aref (fref (mref st-clt) pc-sons) ic) zgmid
                      (ptr (aref lod l)) nofc))
    (cal-bndbox-cog st-clt zgmid lod nofc))

(def (create-cluster nmbr ndpth nstrt nsize ndim nson)
    (fn (ptr cluster) int int int int int int)
  (decl st-clt (ptr cluster))
  (= st-clt (cast (ptr cluster) (csym::malloc (sizeof cluster)))) (= nmbr (+ nmbr 1))
  (= (fref (mref st-clt) nstrt) nstrt) (= (fref (mref st-clt) nsize) nsize)
  (= (fref (mref st-clt) ndim) ndim) (= (fref (mref st-clt) nnson) nson)
  (= (fref (mref st-clt) nmbr) nmbr) (= (fref (mref st-clt) ndpth) ndpth)
  (= (fref (mref st-clt) pc-sons)
     (cast (ptr (ptr cluster)) (csym::malloc (* nson (sizeof (ptr cluster))))))
  (return st-clt))

(def (free-st-clt st-clt) (fn void (ptr cluster))
  (decl ic int)
  (def nnson int (fref (mref st-clt) nnson))
  (for ((= ic 0) (< ic nnson) (inc ic))
    (free-st-clt (aref (fref (mref st-clt) pc-sons) ic)))
  (csym::free (fref (mref st-clt) bmin))
  (csym::free (fref (mref st-clt) bmax))
  (csym::free (fref (mref st-clt) pc-sons)))

(def (checkClusterTree f st-clt) (fn void (ptr FILE) (ptr cluster))
  (if (< (fref (mref st-clt) ndpth) 11)
      (begin
        (csym::fprintf f "%d %d %d %d %lf~%"
                       (fref (mref st-clt) nstrt)
                       (fref (mref st-clt) nsize)
                       (fref (mref st-clt) ndpth)
                       (fref (mref st-clt) nnson)
                       (fref (mref st-clt) zwdth))))
  (if (== (fref (mref st-clt) nnson) 0)
      (begin (return))
      (if (== (fref (mref st-clt) nnson) 1)
          (begin (checkClusterTree f (aref (fref (mref st-clt) pc-sons) 0)))
          (begin (checkClusterTree f (aref (fref (mref st-clt) pc-sons) 0))
                 (checkClusterTree f (aref (fref (mref st-clt) pc-sons) 1))))))

(def (create-ctree-ssgeom st-clt zgmid param lod ndpth ndscd nsrt nd md ndim nclst)
    (fn (ptr cluster) (ptr cluster) (ptr (array double 3)) (array double) (ptr int)
        int int int int int int int)
  (decl id int)
  (decl il int)
  (decl nson int)
  (def minsz double (aref param 21))
  (def zcoef double (aref param 31))
  (decl zlmin (array double ndim))
  (decl zlmax (array double ndim))
  (= ndpth (+ ndpth 1))
  (if (<= nd minsz)
      (begin 
        (= nson 0)
        (= st-clt (create-cluster nclst ndpth nsrt nd ndim nson))
        (if (> ndpth depth-max)
            (begin
              (= depth-max ndpth)))
        (inc count-node))
      (begin
        (for ((= id 0) (< id ndim) (inc id))
          (= (aref zlmin id) (aref (aref zgmid (aref lod 0)) id))
          (= (aref zlmax id) (aref zlmin id))
          (for ((= il 1) (< il nd) (inc il))
            (def zg double (aref (aref zgmid (aref lod il)) id))
            (if (< zg (aref zlmin id)) (begin (= (aref zlmin id) zg))
                (if (< (aref zlmax id) zg) (begin (= (aref zlmax id) zg))))))
        (def zdiff double (- (aref zlmax 0) (aref zlmin 0)))
        (def ncut int 0)
        (for ((= id 0) (< id ndim) (inc id))
          (def zidiff double (- (aref zlmax id) (aref zlmin id)))
          (if (> zidiff (* zcoef zdiff))
              (begin
                (= zdiff zidiff)
                (= ncut id))))
        (def zlmid double (* 0.5 (+ (aref zlmax ncut) (aref zlmin ncut))))
        (def nl int 0) (def nr int (- nd 1))
        (while (< nl nr)
          (while (and (< nl nd)
                      (<= (aref (aref zgmid (aref lod nl)) ncut) zlmid))
            (= nl (+ nl 1)))
          (while (and (>= nr 0)
                      (> (aref (aref zgmid (aref lod nr)) ncut) zlmid))
            (= nr (- nr 1)))
          (if (< nl nr)
              (begin 
                (def nh int (aref lod nl)) (= (aref lod nl) (aref lod nr))
                (= (aref lod nr) nh))))
        (if (or (== nl nd) (== nl 0))
            (begin
              (= nson 1)
              (= st-clt (create-cluster nclst ndpth nsrt nd ndim nson))
              (if (> ndpth depth-max)
                  (begin
                    (= depth-max ndpth)))
              (inc count-node)
              (= (aref (fref (mref st-clt) pc-sons) 0)
                 (create-ctree-ssgeom (aref (fref (mref st-clt) pc-sons) 0) zgmid param
                                      lod ndpth ndscd nsrt nd md ndim nclst)))
            (begin 
              (= nson 2)
              (= st-clt (create-cluster nclst ndpth nsrt nd ndim nson))
              (if (> ndpth depth-max)
                  (begin
                    (= depth-max ndpth)))
              (inc count-node)
              (def nsrt1 int nsrt)
              (def nd1 int nl)
              (= (aref (fref (mref st-clt) pc-sons) 0)
                 (create-ctree-ssgeom (aref (fref (mref st-clt) pc-sons) 0) zgmid param
                                      lod ndpth ndscd nsrt1 nd1 md ndim nclst))
              (= nsrt1 (+ nsrt nl))
              (= nd1 (- nd nl))
              (= (aref (fref (mref st-clt) pc-sons) 1)
                 (create-ctree-ssgeom (aref (fref (mref st-clt) pc-sons) 1) zgmid param
                                      (ptr (aref lod nl)) ndpth ndscd nsrt1 nd1 md ndim nclst))))))
  (= (fref (mref st-clt) ndscd) nd)
  (return st-clt))

(def (csym::get-wall-time) (fn double void)
  (decl time (struct timeval))
  (if (gettimeofday (ptr time) NULL) (begin (return 0)))
  (return (+ (cast double (fref time tv-sec)) (* (cast double (fref time tv-usec)) 0.1))))

(def (csym::get-cpu-time) (fn double void)
  (return (/ (cast double (csym::clock)) CLOCKS-PER-SEC)))