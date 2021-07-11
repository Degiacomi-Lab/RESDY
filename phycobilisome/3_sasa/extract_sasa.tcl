


proc measure_sasa {chain resid fout proberad} {

    set lys_sasa [measure sasa $proberad [atomselect top all] -restrict [atomselect top "chain 0 and resid $resid and sidechain"]]
    set lineout [append $chain " " $resid " " $lys_sasa "\n"]
    puts -nonewline $fout $chain
    puts -nonewline $fout "\t"
    puts -nonewline $fout $resid
    puts -nonewline $fout "\t"
    puts $fout $lys_sasa

}

#S. Sircar,  Basic  Research  Needs  for  Design  of  Adsorptive  Gas  Separation Processes. Ind. Eng. Chem. Res.2006,45(16), 5435-5448
set proberad 1.72

set fout [open "sasa.dat" "w"]


# iterate over selection PDB files (for memory management purposes)
set files [glob ../selections/*.pdb]

foreach f $files {

    mol new $f
    set split1 [split $f /]
    set split2 [split [lindex $split1 2] "."]    
    set chain [lindex $split2 0]

    puts $chain
    set lys_resid [[atomselect top "name CA and resname LYS and chain 0"] get resid]
    foreach resid $lys_resid {
        measure_sasa $chain $resid $fout $proberad
    }

    mol delete top
}

close $fout

