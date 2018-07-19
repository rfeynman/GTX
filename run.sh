#!/usr/bin/tcsh
echo input $1
set p = $1
#echo $p
@ von = $2
if ( $# == 3 ) then
	@ bis = $3
else
	@ bis = $von
endif

#echo $p $von $bis $#
@ run = $von

set eins = 1
while ( $run <= $bis )
	echo gpt  -o outgpt.out0$run $p section=$run
	gpt  -o outgpt.out0$run $p section=$run
	if ( $status  ) then
        	echo "aborted..."
        	exit
	endif
	echo "completed $run ..."
	@ run = $run + $eins
	#echo $run $bis
end




#~/GTX/gplot outgpt.out0[0123]
#xmgrace -geometry 780x630+600 -fixed 510 380 -noask env.agr &
