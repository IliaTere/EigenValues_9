for ((n=1; n<=12;n++)) ; do for ((s=3;s<=3;s+=1)) ; do echo $'\n' "n=$n s=$s --------------" ; ./a.out $n 0 1e-14 $s ; ./a.out $n 0 1e-15 $s; ./a.out $n 0 1e-13 $s ; done ; done
