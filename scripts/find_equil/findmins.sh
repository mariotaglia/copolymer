for i in curvature* ;
	do cd $i
	for j in ph* ;
	    do cd $j 
		python3 ~/develop/copolymer/scripts/find_equil/pick_eq.py
                cd equil
	        ~/develop/copolymer/assembly
	    cd ..
	  cd ..
        done
cd ..
done
