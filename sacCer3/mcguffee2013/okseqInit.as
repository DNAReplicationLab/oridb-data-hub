table okseqInit
"Okazaki-fragment initiation zones (McGuffee et al. 2013, re-analysed on sacCer3)"
    (
    string chrom;          "Reference sequence chromosome"
    uint   chromStart;     "Start of the initiation zone"
    uint   chromEnd;       "End of the initiation zone"
    string name;           "Zone identifier"
    uint   score;          "1000 x origin efficiency (OEM)"
    char[1] strand;        "Not applicable"
    uint   thickStart;     "Derivative peak position"
    uint   thickEnd;       "Derivative peak position + 1"
    float  efficiency_OEM; "Origin efficiency metric, raw (not deconvolved)"
    float  delta_FD;       "Change in fork directionality across the zone"
    float  deriv_per_kb;   "Peak derivative of fork directionality, per kb"
    uint   width_bp;       "Zone width in bp"
    )
