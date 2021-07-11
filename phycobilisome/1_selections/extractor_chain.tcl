# get pdbs of chains as well as adjacent chains


proc extractor {chain} {

   set chlist {a b c d e f g h i j k l m n o p q r s t u v w 1 2 3 4 5 6 7 8 9}

   set ext ".pdb"
   set fout $chain$ext
   set ch_neigh [[atomselect top "protein and (chain $chain or (within 5 of (protein chain $chain)))"] get chain]
   set ch_neigh_unique [lsort -unique $ch_neigh]

   #problem: PDB format only allows single letter code for chains. Rename them all!

   set cnt 0
   set mychains {}
   foreach c $ch_neigh_unique {

       set a [atomselect top "protein and chain $c"]
       if {$c==$chain} {
           $a set chain 0
           lappend mychains 0
       } else {
           set thischain [lindex $chlist $cnt]
           $a set chain $thischain
           incr cnt
           lappend mychains $thischain
       }
   }

   # save all relevant chains in file of desired name
   [atomselect top "(protein or (resname MEN and not name CE2)) and chain $mychains"] writepdb $fout

}

###################################


# get chains containing at least one lysine
mol new ../6kgx-assembly1.cif
set lys_chain [[atomselect top "name CA and resname LYS"] get chain]
set lys_chain_unique [lsort -unique $lys_chain]
mol delete top

set lys_chain_unique {LW NX}

# iterate over all these chains
foreach i $lys_chain_unique {
   mol new ../6kgx-assembly1.cif
   puts $i

   set code [catch {
        extractor $i
   } result]

   if {$code > 0} {
       puts "failed!"
   } 


   mol delete top
}

