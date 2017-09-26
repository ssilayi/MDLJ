#!bin/bash

#double ti[11] = {2.0,1.5,1.0,0.8,0.75,0.7,0.65,0.6,0.5,0.3,0.1};
#double ri[2] = {1.0, 0.7};
reset
set term pdf
set output 'rdf_1.pdf'

set key outside
set xlabel "distance r"
set ylabel "g(r)"

plot '1/rdf_0.1.txt' u 1:2 w l t '{/Symbol r} = 1, T = 0.1',  \
     '1/rdf_0.3.txt' u 1:2 w l t '{/Symbol r} = 1, T = 0.3',  \
     '1/rdf_0.5.txt' u 1:2 w l t '{/Symbol r} = 1, T = 0.5',  \
     '1/rdf_0.6.txt' u 1:2 w l t '{/Symbol r} = 1, T = 0.6',  \
     '1/rdf_0.65.txt' u 1:2 w l t '{/Symbol r} = 1, T = 0.65',\
     '1/rdf_0.7.txt' u 1:2 w l t '{/Symbol r} = 1, T = 0.7',  \
     '1/rdf_0.75.txt' u 1:2 w l t '{/Symbol r} = 1, T = 0.75', \
     '1/rdf_0.8.txt' u 1:2 w l t '{/Symbol r} = 1, T = 0.8',  \
     '1/rdf_1.txt' u 1:2 w l t '{/Symbol r} = 1, T = 1.0',  \
     '1/rdf_1.5.txt' u 1:2 w l t '{/Symbol r} = 1, T = 1.5',  \
     '1/rdf_2.txt' u 1:2 w l t '{/Symbol r} = 1, T = 2.0'

    
reset
set term pdf
set output 'rdf_0.7.pdf'

set key outside
set xlabel "distance r"
set ylabel "g(r)"

plot '0.7/rdf_0.1.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 0.1',  \
     '0.7/rdf_0.3.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 0.3',  \
     '0.7/rdf_0.5.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 0.5',  \
     '0.7/rdf_0.6.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 0.6',  \
     '0.7/rdf_0.65.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 0.65', \
     '0.7/rdf_0.7.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 0.7',  \
     '0.7/rdf_0.75.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 0.75', \
     '0.7/rdf_0.8.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 0.8',  \
     '0.7/rdf_1.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 1.0',  \
     '0.7/rdf_1.5.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 1.5',  \
     '0.7/rdf_2.txt' u 1:2 w l t '{/Symbol r} = 0.7, T = 2.0'
