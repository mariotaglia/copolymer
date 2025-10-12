listapH='7.'


for i in `cat lista` ; do
        for j in $listapH ; do
    echo $i $j
    python3 ~/develop/copolymer/scripts/auto/PA.py $i $j
done
done

