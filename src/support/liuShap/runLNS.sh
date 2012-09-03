mkdir outputFiles
mkdir images
mkdir images/sigma
mkdir images/FJ

rm outputFiles/*
rm images/*/*.png

./LNS > analytic.dat
./Tango -i analytic.dat > tango.out
python plot.py

mv images/*sigma.png images/sigma/
mv images/*.png images/FJ/
