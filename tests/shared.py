"""Shared constants for crisviper tests."""

# Standard CARLIN 332bp reference sequence
CARLIN_REF = (
    "TATGTGTGGGAGGGCTAAGAGG"  # Primer5 (23bp)
    "CCGCC"                    # Prefix (5bp)
    "GACTGCACGACAGTCGA"        # Target1 (13+7 cutsite)
    "CGATGGAG"                 # Linker (7bp)
    "TCGACACGACTCGCGCA"        # Target2
    "TACGATGG"                 # Linker
    "AGTCGACTACAGTCGCTA"       # Target3
    "CGACGATG"                 # Linker
    "GAGTCGCGAGCGCTATG"        # Target4
    "AGCGACTA"                 # Linker
    "TGGAGTCGATACGATACG"       # Target5
    "CGCACGCT"                 # Linker
    "ATGGAGTCGAGAGCGCGC"       # Target6
    "TCGTCAAC"                 # Linker
    "GATGGAGTCGCGACTGTA"       # Target7
    "CGCACTCG"                 # Linker
    "CGATGGAGTCGATAGTAT"       # Target8
    "GCGTACAC"                 # Linker
    "GCGATGGAGTCGACTGCA"       # Target9
    "CGACAGTC"                 # Linker
    "GACTATGGAGTCGATACGTAGC"   # Target10
    "ACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"  # Postfix(8)+Primer3(33)
)
