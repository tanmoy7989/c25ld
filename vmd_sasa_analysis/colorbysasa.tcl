### USAGE:
## For CG:-  vmd -dispdev win -e colorbysasa.tcl -args <lammpstraj> 0 25 0 20000 10 0 
## For AA:-  vmd -dispdev win -e colorbysasa.tcl -args <lammpstraj> 5100 5125 0 20000 100 3


# Coloring routine
proc colorbysasa {polymer_sel serial_start serial_stop} {
	set mon_rad 2.0934
	set w_rad 1.4
	$polymer_sel set radius $mon_rad
	for {set j $serial_start} {$j < $serial_stop} {incr j} {	
		set monomer_sel [atomselect top "index $j"]				
		set sasa [measure sasa $w_rad $polymer_sel -restrict $monomer_sel]
		$monomer_sel set user $sasa
		$monomer_sel delete
	}
	mol modcolor 0 [molinfo top] User
	#mol colupdate 0 [molinfo top] 1
	mol scaleminmax [molinfo top] 0 auto
}

### MAIN ####

# Parse user input
set trajfile   	 	[lindex $argv 0]
set serial_start 	[lindex $argv 1]
set serial_end   	[lindex $argv 2]
set start      		[lindex $argv 3]
set stop       		[lindex $argv 4]
set freq       		[lindex $argv 5]
set mon_type		[lindex $argv 6]

# Read in files
set maxframes  [expr int(($stop-$start)/$freq)]	
mol new $trajfile first $start last $stop step $freq waitfor $maxframes

set mon_rad 2.0934
set w_rad 1.4
mol modstyle 0 top VDW
color Display Background white
display projection orthographic

set nframes [molinfo top get numframes]
set polymer_sel [atomselect top "type $mon_type"]
$polymer_sel set radius $mon_rad

# SASA calculation loop
for {set i 0} {$i <= $nframes} {incr i $freq} {
	molinfo top set frame $i
	puts "Processing frame $i..."
	colorbysasa $polymer_sel $serial_start $serial_end
}

