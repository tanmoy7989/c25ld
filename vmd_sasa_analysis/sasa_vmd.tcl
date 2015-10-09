# Parse user input
set trajfile   	 	[lindex $argv 0]
set outprefix 		[lindex $argv 1]
set serial_start 	[lindex $argv 2]
set serial_end   	[lindex $argv 3]
set start      		[lindex $argv 4]
set stop       		[lindex $argv 5]
set freq       		[lindex $argv 6]
set mon_type		[lindex $argv 7]

# Read in files
set maxframes  [expr int(($stop-$start)/$freq)]	
mol new $trajfile first $start last $stop step $freq waitfor $maxframes
set output [open ${outprefix}.dat w]

set mon_rad 2.0934
set w_rad 1.4

puts " starting vmd for $outprefix..."
set nframes [molinfo top get numframes]
set polymer_sel [atomselect top "type $mon_type"]
$polymer_sel set radius $mon_rad

# SASA calculation loop
for {set i 0} {$i <= $nframes} {incr i $freq} {
	molinfo top set frame $i
	for {set j $serial_start} {$j < $serial_end} {incr j} {	
		set monomer_sel [atomselect top "index $j"]				
		set sasa [measure sasa $w_rad $polymer_sel -restrict $monomer_sel]
		puts -nonewline $output "$sasa\t"
	}
	puts $output "\n"
}
close $output
quit
