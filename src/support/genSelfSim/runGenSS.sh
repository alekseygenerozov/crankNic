function run_single
{
	sed s/inner_bndry_value\ =\ $1/inner_bndry_value\ =\ $2/ params.in > tmp
	mv tmp params.in
	python makeSelfSim.py > selfSim.init
	./Tango selfSim.init
	cp outputFiles/T003.dat savedOutputs/chi_$3_$4.dat
}

function run_diff_case
{

	echo ""
	echo ""
	echo "====================================="
	echo "Running for Test Case $5"

	rm outputFiles/*

	sed s/np\ =\ $1/np\ =\ $2/ params.in > tmp
	sed s/nd\ =\ $3/nd\ =\ $4/ tmp > tmp2
	mv tmp2 params.in
	rm tmp
	
	run_single 1.0 0.0 p0 $5
	run_single 0.0 0.2 p2 $5
	run_single 0.2 0.5 p5 $5
	run_single 0.5 1.0 1p $5

	cp params.out savedOutputs/params_$5.out
}

rm images/*.png
rm savedOutputs/*
run_diff_case -.8 -1.2 .3 .4 1
run_diff_case -1.2 -.8 .4 .3 2
python plotGenSelfSim.py
