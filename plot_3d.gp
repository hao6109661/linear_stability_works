reset
# Set the size of the plot
#set terminal pngcairo size 600,600
set term qt
    
# Set the output file 
#set output sprintf("plot2.png")
    
file1 = sprintf("RESLT/trace.dat")

# Got the number of I values in the general file 
stats file1 nooutput
num_lines = STATS_records

filename = sprintf("data.dat")
set print filename 
    
# Loop over all the files 
do for [i=0:num_lines-1] {

# Set file handle
file2 = sprintf("RESLT/eigenvalues%d.dat", i)

#eigen = system(sprintf("awk '{if(NR==1 || $1 > max) max=$1} END {print max}' \"%s\"", file2))
#print eigen

# The number of the positive eigenvalues 
count = system(sprintf("awk 'BEGIN {count=0} $1 > 1.0e-3 {count++} END {print count}' \"%s\"", file2))
#print count


I = system(sprintf("awk 'NR==%d+1 {print $1}' %s",i, file1))
#print t

theta_eq = system(sprintf("awk 'NR==%d+1 {print $3}' %s",i, file1))

print sprintf("%s %s %s", I, theta_eq, count)


}
unset print

#set style line 1 lc rgb "blue" lt 1 lw 2

#set log x

set label 1 sprintf("q=0.4,alpha=22.5") at screen 0.99,0.03 right front font ",15"

# Define label size and offset
label_size = 15
label_offset_x = 1
label_offset_y = 0.7

# Set the plot range
set xrange [*:*]
set yrange [*:*]
set zrange [*:5]

# Set labels for the axes
set xlabel "I"
set ylabel "max-eigenvalue"

set view 0, 0       
set xlabel "X"
set ylabel "Y"
set zlabel "Z"

set grid               
set ticslevel 0        
set pointsize 1.5      
set style data points  

splot "data.dat" using 1:2:3 with points pointtype 7 pointsize 0.5 linecolor rgb "blue"




# Save the final multiplot layout to the output file
set output

#system(sprintf("rm data.dat"))
