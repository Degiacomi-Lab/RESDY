cat $1 | grep -v " H$" | sed 's/BLEU/ LEU/g' | sed 's/BSER/ SER/g'  | sed 's/BASN/ ASN/g'  | sed 's/BHIS/ HIS/g'|\
    sed 's/BLYS/ LYS/g' | sed 's/BGLU/ GLU/g'  | sed 's/BPRO/ PRO/g' | sed 's/BILE/ ILE/g'|\
    sed 's/BCYS/ CYS/g' | sed 's/BTYR/ TYR/g' | sed 's/BTHR/ THR/g' | sed 's/BASP/ ASP/g' |\
    sed 's/BMET/ MET/g' | sed 's/BGLN/ GLN/g' | sed 's/BGLY/ GLY/g' | sed 's/BTRP/ TRP/g' |\
    sed 's/BARG/ ARG/g' | sed 's/BVAL/ VAL/g' | sed 's/BALA/ ALA/g' | sed 's/BPHE/ PHE/g'|\
    grep -vE "[AC]ALA|[AC]PHE|[AC]HIS|[AC]LEU|[AC]SER|[AC]ASN|[AC]LYS|[AC]GLU|[AC]PRO|[AC]ILE|[AC]CYS|[AC]TYR|[AC]THR|[AC]ASP|[AC]MET|[AC]GLN|[AC]GLY|[AC]TRP|[AC]ARG|[AC]VAL" > tmp