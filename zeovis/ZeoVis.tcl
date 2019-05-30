# VMD general scene settings
mol default color element
# turn on lights 0 and 1
light 0 on
light 1 on
light 2 off
light 3 off
# position the stage and axes
axes location off
stage location off

# Background to white
color Display Background white

# Start of OS-dependent functions
proc remove_file {fileName} {
    exec rm $fileName
}

proc copy_file {oldName newName} {
    exec cp $oldName $newName
}

proc run_network_program {list_of_args} {
    eval exec network $list_of_args
}
# End of OS-dependent functions

proc readFileName {} {
    global filename
    puts "\nPlease enter the filename of your input file\n"
    gets stdin filename
    return ""
}

#proc chooseAccuracy {} {
#  global acculevel
#  puts"Do you want to proceed with default setting or higher accuracy?\[Yes/No\]"
#  gets stdin acculevel
#}

proc readAccuFlag {} {
    global accuflag
    puts "Please enter the desired accuracy for Voronoi decomposition\n"
    puts "Available options:\n"
    puts "None(means running zeo++ without -ha flag)"
    puts "HI MED LOW DEF"
    puts "S4 S10 S20 S30 S40 S50 S100 S500 S1000 S10000"
    puts "FCC ACC AQC DDH TIH ICH ICC RIH\n"
    gets stdin accuflag
    return ""
}

proc checkInputFile {} {
    global filename
    global ending
    if {[string match "*.cuc" $filename]} {
	set ending ".cuc"
    } elseif {[string match "*.cssr" $filename]} {
	set ending ".cssr"
    } elseif {[string match "*.cif" $filename]} {
	set ending ".cif"
    } else {
	puts "ERROR: INPUT FILE MUST BE OF TYPE .cuc or .cssr \nExiting"
	exit
    }
}


readFileName
readAccuFlag
checkInputFile
set inputFile "ZeoVisInput$ending"
copy_file $filename $inputFile
puts "CREATING zeovis INPUT FILE"

#run_network_program [list "-r" "example.rad" "-zvis" "ZeoVisInput.zvis" $inputFile]
if {[string equal $accuflag "None"]} {
  run_network_program [list "-r" "-zvis" "ZeoVisInput.zvis" $inputFile]
} else {
  run_network_program [list "-ha" $accuflag "-r" "-zvis" "ZeoVisInput.zvis" $inputFile]
}

#run_network_program [list "-r" "-ha" "-zvis" "ZeoVisInput.zvis" $inputFile]
#run_network_program [list "-zvis" "ZeoVisInput.zvis" $inputFile]
puts "INPUT FILE SUCCESSFULLY CREATED"

source "ZeoVisInput.zvis"
remove_file "ZeoVisInput.zvis"

source "ZeoVisCommands.tcl"
set num_segments 0
set num_features 0
set num_cages    0
initialize
puts "ENJOY ZEOVIS!"
