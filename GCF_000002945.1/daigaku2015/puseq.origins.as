table puseqOrigins
"Daigaku 2015 Pu-seq origin calls on ASM294v2"
    (
    string  chrom;       "Reference sequence"
    uint    chromStart;  "Start, 0-based"
    uint    chromEnd;    "End, exclusive"
    string  name;        "Call identifier, chromosome and position"
    uint    score;       "Efficiency scaled to 0-1000 over the published range"
    char[1] strand;      "Not applicable"
    string  efficiency;  "Published efficiency, percentage of cells"
    )
