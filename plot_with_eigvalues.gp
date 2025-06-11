reset
# Set the size of the plot
#set terminal pngcairo size 600,600
set term qt
    
# Set the output file 
#set output sprintf("plot2.png")
    
file1 = sprintf("RESLT_q_0.4_alpha_22.5_initial_-4.80/trace.dat")

# Got the number of I values in the general file 
stats file1 nooutput
num_lines = STATS_records

filename = sprintf("data1.dat")
set print filename 
    
# Loop over all the files 
do for [i=0:num_lines-1] {

# Set file handle
file2 = sprintf("RESLT_q_0.4_alpha_22.5_initial_-4.80/eigenvalues%d.dat", i)

eigen = system(sprintf("awk '{if(NR==1 || $1 > max) max=$1} END {print max}' \"%s\"", file2))
#print eigen

I = system(sprintf("awk 'NR==%d+1 {print $1}' %s",i, file1))
#print t

theta_eq = system(sprintf("awk 'NR==%d+1 {print $3}' %s",i, file1))

print sprintf("%s %s %s", I, theta_eq, eigen)


}
unset print

#set style line 1 lc rgb "blue" lt 1 lw 2

set label 1 sprintf("q=0.3,alpha=22.5") at screen 0.99,0.03 right front font ",15"

# Define label size and offset
label_size = 15
label_offset_x = 1
label_offset_y = 0.7

# Set the plot range
set xrange [*:*]
set yrange [*:*]

# Set labels for the axes
set xlabel "I"
set ylabel "max-eigenvalue"

plot "data1.dat" using 1:2 with p lc rgb "red" title "", \
     "data1.dat" using 1:2:3 with labels offset label_offset_x,label_offset_y font ",label_size" tc rgb "black" title "", \
     "../plots/RESLT_q_0.30_alpha_0.125pi_element_20/elastic_beam_I_theta_q_0.30_alpha_0.125pi_initial_-4.80.dat" u 1:3 w l


# Save the final multiplot layout to the output file
set output

#system(sprintf("rm data.dat"))
