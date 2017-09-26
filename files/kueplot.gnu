#!bin/bash


#1-4   T_{set}   T_{sys}  T_{ave}   	T_{var}
#5-7   K/N   	K_{ave}/N	K_{var}
#8-10  U/N   	U_{ave}/N 	U_{var}
#11-13 E/N   	E_{ave}/N 	E_{var}
#14-16 P     	P_{ave}   	P_{var}
# tfile << T << "\t" << K_ave[cc-1] / double (N) << "\t" << U_ave[cc-1] / double (N) << "\t" 
#                  << E_ave[cc-1]/ double (N) << "\t" << p_ave[cc-1]<< "\t" << t_ave[cc-1]
#                  << "\t U_corr : " << ucor <<"\t P_corr : " << pcor <<"\t runtime:" << sf-si << "\n";

reset
set term pdf
set key outside
set output "kvst.pdf"
set title "K/N vs. T"
set xlabel "T"
set ylabel "K/N"
plot "1/tracking.txt" u 1:6 w lp pt 7 t '{/Symbol r} = 1', \
     "0.7/tracking.txt" u 1:6 w lp pt 8 t '{/Symbol r} = 0.7'

reset            
set term pdf
set key outside
set output "uvst.pdf"
set title "U/N vs. T"
set xlabel "T"
set ylabel "U/N"
plot "1/tracking.txt" u 1:9 w lp  pt 7 t '{/Symbol r} = 1', \
     "0.7/tracking.txt" u 1:9 w lp pt 8 t '{/Symbol r} = 0.7'


reset
set term pdf
set key outside
set output "evst.pdf"
set title "E/N vs. T"
set xlabel "T"
set ylabel "E/N"
plot "1/tracking.txt" u 1:12 w lp pt 7t '{/Symbol r} = 1', \
     "0.7/tracking.txt" u 1:12 w lp pt 8 t '{/Symbol r} = 0.7'

reset
set term pdf
set key outside
set output "pvst.pdf"
set title "P vs. T"
set xlabel "T"
set ylabel "P"
plot "1/tracking.txt" u 1:15 w lp pt 7 t '{/Symbol r} = 1', \
     "0.7/tracking.txt" u 1:15 w lp pt 8 t '{/Symbol r} = 0.7'

reset
set term pdf
set key outside
set output "tvst.pdf"
set title "T vs. T"
set xlabel "T"
set ylabel "P"
plot "1/tracking.txt" u 1:3 w lp pt 7 t '{/Symbol r} = 1', \
     "0.7/tracking.txt" u 1:3 w lp pt 8 t '{/Symbol r} = 0.7'
