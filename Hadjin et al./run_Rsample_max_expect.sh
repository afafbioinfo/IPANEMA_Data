for i in "5domainof16SrRNAEcoli" "5domainof16SrRNAHvolcanii" "5domainof23SrRNAEcoli"  "5SrRNAEcoli"  "cyclicdiGMPriboswitchVcholerae" "AdenineriboswitchVvulnificus" "GroupIIntronTthermophila"  "HepatitisCvirusIRESdomain"  "HIV15primepseudoknotdomain"  "LysineriboswitchTmaritime"  "P546domainbI3groupIintron"  "PreQ1riboswitchBsubtilis"  "RNasePBsubtilis"  "SAMIriboswitchTtengcongensis"  "tRNAaspyeast"  "MBoxriboswitchBsubtilis" "FluorideriboswitchPsyringae" "SARScoronaviruspseudoknot"  "TPPriboswitchEcoli" "GroupIIintronOiheyensis" "GroupIintronAzoarcussp" "SignalrecognitonparticleRNAhuman" "Telomerasepseudoknothuman" "tRNApheEcoli"

do
   SECONDS=0
        echo ${i}
        exe/Rsample fasta_files/$i.fa  Shape_files/$i"1M71Shape.txt"  24newpart/$i.pfs 
        exe/MaxExpect 24newpart/$i.pfs  max_expect_ctnew/$i.ct -s 1
        echo $SECONDS


done
