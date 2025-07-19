set terminal x11 
set title 'Lorenz Attractor - Animation'
set xlabel 'X'; set ylabel 'Y'; set zlabel 'Z'

do for [i=1:100000:10] {
    splot 'solution.dat' every ::1::i using 2:3:4 with lines notitle
    pause 0.01
}
