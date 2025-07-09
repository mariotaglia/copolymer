for i in curvature* ;
	do cd $i
	for j in ph* ;
	    do cd $j 
		python3 ~/develop/copolymer/scripts/find_equil/pick_eq.py
	        ~/develop/copolymer/assembly
	    cd ..
        done
cd ..
done
