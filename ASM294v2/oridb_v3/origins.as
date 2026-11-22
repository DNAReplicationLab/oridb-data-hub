table oridbOrigins
"OriDB 3.0.0-trial.1 replication origins"
    (
    string  chrom;        "Reference sequence"
    uint    chromStart;   "Start, 0-based"
    uint    chromEnd;     "End, exclusive"
    string  name;         "OriDB accession"
    uint    score;        "Shaded by tier: Confirmed 1000 to Not assessed 200"
    char[1] strand;       "Not applicable to an origin"
    uint    thickStart;   "Same as chromStart"
    uint    thickEnd;     "Same as chromEnd"
    uint    reserved;     "Tier colour, as itemRgb"
    string  tier;         "Evidence tier"
    string  arsNegative;  "Inside a fragment that failed an ARS assay"
    uint    nSupport;     "Datasets that support this site"
    string  originName;   "Historical OriDB name, where one exists"
    string  okseqOem;     "Ok-seq firing efficiency (OEM), or '.' if unresolved"
    string  replTime;     "Median replication time in minutes, or '.' if unresolved"
    )
