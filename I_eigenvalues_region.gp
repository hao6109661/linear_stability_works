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

# Initialize an empty plot command string
cmd = ""
    
# Loop over all the files 
do for [i=1680:1700] {

# Set file handle
file2 = sprintf("RESLT/eigenvalues%d.dat", i)

I = system(sprintf("awk 'NR==%d+1 {print $1}' %s",i, file1))

    # For the first entry, do not add a comma
    if (i == 1680) {
        cmd = sprintf('"%s" using (%s):1 with linespoint pointsize 1 title ""', file2, I)
    } else {
        # For subsequent entries, add a comma before appending
        cmd = sprintf('%s, "%s" using (%s):1 with linespoint pointsize 1 title ""', cmd, file2, I)
    }

}

# Set x and y axis range to auto
set xrange [*:*]
set yrange [-0.1:0.1]

# Execute the full plot command
eval("plot " . cmd)



# Save the final multiplot layout to the output file
set output

