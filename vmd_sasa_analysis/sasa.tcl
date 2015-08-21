# Parse user input
set trajfile   	 	[lindex $argv 0]
set outprefix 		[lindex $argv 1]
set serial_start 	[lindex $argv 2]
set serial_end   	[lindex $argv 3]
set start      		[lindex $argv 4]
set stop       		[lindex $argv 5]
set freq       		[lindex $argv 6]

# Read in files
set maxframes  [expr int(($stop-$start)/$freq)]	
mol new $trajfile first $start last $stop step $freq waitfor $maxframes

set nframes [molinfo top get numframes]
set output [open ${outprefix}.dat w]
set polymer_sel [atomselect top all]

# SASA calculation loop
for {set i 0} {[expr $i*$freq] <= $nframes} {incr i} {
	molinfo top set frame $i
	for {set j $serial_start} {$j < $serial_end} {incr j} {	
		set monomer_sel [atomselect top "index $j"]
		set sasa [measure sasa 1.4 $polymer_sel -restrict $monomer_sel]
		puts -nonewline $output "$sasa\t"
	}
	puts $output "\n"
}
close $output
quit
