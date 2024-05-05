## set the dummy position for the extra atoms - (xd,yd,zd)
set xd 0.0
set yd 0.0
set zd 0.0
set env(LAMMPSDUMMYPOS) {$xd,$yd,$zd}

## set the maximum possible number of atoms, Nmax
## - this needs to be edited for specific trajectories
set Nmax 20000
set env(LAMMPSMAXATOMS) $Nmax

## remap fields in lammpstrj to velocities with LAMMPS plugin
set env(LAMMPSREMAPFIELDS) {vx=c_id_track,vy=c_type_track}

# Initial Setup
color Display {Background} white
color scale method viridis
color scale min 0.0
color scale midpoint 0.5
color scale max 1.0
display projection Orthographic

# Define a switch for the coloring methodology
#  false - color based on instantaneous topology
#  true - color based on final topology
set final_topo_color false

# Define Radii
set DNA_rad 13.0
set a_special 2.0
set ori_rad [expr $a_special*$DNA_rad]
set ter_rad [expr $a_special*$DNA_rad]
set fork_rad [expr $a_special*$DNA_rad]
set hinge_rad [expr 1.5*$DNA_rad]
set anchor_rad [expr 1.5*$DNA_rad]
set ribo_rad 70.0
set bdry_rad 32.5

# set the file to be loaded
set molecule_file test_trajectory.lammpstrj

# load the new molecule
mol new $molecule_file  waitfor all type lammpstrj step 1
set mol_id [molinfo top]
set nframes [molinfo $mol_id get numframes]

# jump to the end
animate goto end

# prepare atom selections
set bdry [atomselect $mol_id "vy==1 and z<0" frame now]
set ribo [atomselect $mol_id "vy==2" frame now]
set DNA [atomselect $mol_id "vy>2" frame now]
set mono [atomselect $mol_id "vy==3" frame now]
set ori [atomselect $mol_id "vy==4" frame now]
set ter [atomselect $mol_id "vy==5" frame now]
set fork [atomselect $mol_id "vy==6" frame now]
set anchor [atomselect $mol_id "vy==7" frame now]
set hinge [atomselect $mol_id "vy==8" frame now]

# Center and zoom the loaded molecule
mol modstyle 0 $mol_id Points
# mol modcolor 0 0 ColorID 1
mol drawframes $mol_id 0
mol representation Points 10.0
mol delrep 0 $mol_id

# create a counter for the number of representations
set rep_count 0

# create representation for DNA
mol addrep $mol_id
set DNA_rep $rep_count
set rep_count [expr $rep_count + 1]
mol modstyle $DNA_rep $mol_id VDW $DNA_rad 12.0
mol modmaterial $DNA_rep $mol_id AOChalky
mol modselect $DNA_rep $mol_id [$DNA text]
mol scaleminmax $mol_id $DNA_rep 1 [$DNA num]
mol selupdate $DNA_rep $mol_id on
for {set i 0} {$i<$nframes} {incr i} {
    animate goto $i
    #$DNA frame $i
    $DNA update
    puts "Setting User data for frame [$DNA frame] with [$DNA num] monomers ..."
    if {$final_topo_color == true} {
	set idx_temp [$DNA get vx]
    } else {
	set n_DNA [$DNA num]
	set idx_temp [lmap x [$DNA get vx] {expr $x/$n_DNA}]
    }
    $DNA set user $idx_temp
}
mol modcolor $DNA_rep $mol_id User
mol colupdate $DNA_rep $mol_id on

# create representation for oris
mol addrep $mol_id
set ori_rep $rep_count
set rep_count [expr $rep_count + 1]
mol modstyle $ori_rep $mol_id VDW $ori_rad 20.0
mol modmaterial $ori_rep $mol_id AOChalky
mol modselect $ori_rep $mol_id [$ori text]
mol selupdate $ori_rep $mol_id on
mol modcolor $ori_rep $mol_id ColorID 1
mol colupdate $ori_rep $mol_id on

# create representation for ters
mol addrep $mol_id
set ter_rep $rep_count
set rep_count [expr $rep_count + 1]
mol modstyle $ter_rep $mol_id VDW $ter_rad 20.0
mol modmaterial $ter_rep $mol_id AOChalky
mol modselect $ter_rep $mol_id [$ter text]
mol selupdate $ter_rep $mol_id on
mol modcolor $ter_rep $mol_id ColorID 3
mol colupdate $ter_rep $mol_id on

# create representation for forks
mol addrep $mol_id
set fork_rep $rep_count
set rep_count [expr $rep_count + 1]
mol modstyle $fork_rep $mol_id VDW $fork_rad 20.0
mol modmaterial $fork_rep $mol_id AOChalky
mol modselect $fork_rep $mol_id [$fork text]
mol selupdate $fork_rep $mol_id on
mol modcolor $fork_rep $mol_id ColorID 27
mol colupdate $fork_rep $mol_id on

# create representation for anchors
mol addrep $mol_id
set anchor_rep $rep_count
set rep_count [expr $rep_count + 1]
mol modstyle $anchor_rep $mol_id VDW $anchor_rad 12.0
mol modmaterial $anchor_rep $mol_id AOChalky
mol modselect $anchor_rep $mol_id [$anchor text]
mol selupdate $anchor_rep $mol_id on
mol modcolor $anchor_rep $mol_id ColorID 16
mol colupdate $anchor_rep $mol_id on

# create representation for hinges
mol addrep $mol_id
set hinge_rep $rep_count
set rep_count [expr $rep_count + 1]
mol modstyle $hinge_rep $mol_id VDW $hinge_rad 12.0
mol modmaterial $hinge_rep $mol_id AOChalky
mol modselect $hinge_rep $mol_id [$hinge text]
mol selupdate $hinge_rep $mol_id on
mol modcolor $hinge_rep $mol_id ColorID 8
mol colupdate $hinge_rep $mol_id on

# create representation for ribosomes
mol addrep $mol_id
set ribo_rep $rep_count
set rep_count [expr $rep_count + 1]
mol modstyle $ribo_rep $mol_id VDW $ribo_rad 20.0
mol modmaterial $ribo_rep $mol_id AOChalky
mol modselect $ribo_rep $mol_id [$ribo text]
mol selupdate $ribo_rep $mol_id on
mol modcolor $ribo_rep $mol_id ColorID 13
mol colupdate $ribo_rep $mol_id on


# create representation for boundary particles
mol addrep $mol_id
set bdry_rep $rep_count
set rep_count [expr $rep_count + 1]
mol modstyle $bdry_rep $mol_id VDW $bdry_rad 20.0
mol modmaterial $bdry_rep $mol_id Transparent
mol modselect $bdry_rep $mol_id [$bdry text]
mol selupdate $bdry_rep $mol_id on
mol modcolor $bdry_rep $mol_id ColorID 2
mol colupdate $bdry_rep $mol_id on


# Center the view
display resetview
