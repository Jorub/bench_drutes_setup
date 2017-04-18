#!/bin/bash

function run_drutes {

 # rm -rf $1
  #mkdir $1
  #cp -a drutemp/* $1

  cd $1
  
  alpha=$2
  n=$3
  m=`echo "scale=12; 1.0-1.0/$n" | bc`
  ths=$4
  Ks=$5

  #   substitution of parameters into input files for drutes 
  sed -e 's/!alpha/'$alpha'/g' -e 's/!n/'$n'/g' -e 's/!m/'$m'/g' -e 's/!ths/'$ths'/g' -e 's/!Ks/'$Ks'/g' drutes.conf/water.conf/matrixinv.conf > drutes.conf/water.conf/matrix.conf
  bin/drutes > /dev/null
  echo 'process': $1 'done' > objfnc.val
  cd ..
    
    
}


rm -f drutes.vals

#count the number of processes       
let nproc=0
while read l a b c d; do
    if [[  $l == "p"  ]]; then
      let nproc=nproc+1
    fi
  done < pars.in
  
#execute drutes function in parallel
let z=0
while read l a b c d
  do
    if [[  $l == "p"  ]]; then
      let z=z+1
      if [[ $z -lt $nproc ]] ; then
        run_drutes $z $a $b $c $d &
      else
        run_drutes $z $a $b $c $d
      fi
    fi
  done < pars.in
 

#at the end of drutes function file obj.val is created, if all files for each process exist then we have finished 
let z=0    
while [[ $z -lt $nproc ]]
  do
  let z=0
    for i in `seq 1 $nproc`
      do
        FILE="$i/objfnc.val"
        if [ -f $FILE ];
          then
             let z=z+1
        fi
    done    
done

for i in `seq 1 $nproc` ; do
  val=`cat $i/objfnc.val`
  echo $i $val >> drutes.vals
done
