table trepPeaks
"Replication origin calls from the Muller et al. (2014) sync-seq Trep profile (local minima of replication time)"
    (
    string  chrom;        "Reference sequence"
    uint    chromStart;   "Start of the call, 0-based"
    uint    chromEnd;     "End of the call: the band within 1 minute of the trough minimum"
    string  name;         "Call identifier, numbered in genome order"
    uint    score;        "Prominence scaled to 0-1000 for shading"
    char[1] strand;       "Not applicable"
    uint    thickStart;   "Summit: midpoint of the 1 kb bin holding the earliest Trep"
    uint    thickEnd;     "Summit end (thickStart + 1)"
    float   trep;         "Smoothed Trep at the summit, minutes from alpha-factor release"
    float   prominence;   "Depth of the trough below its enclosing saddle, minutes"
    uint    interpolated; "1 if an unfitted bin was interpolated within 1 kb of the summit"
    )
