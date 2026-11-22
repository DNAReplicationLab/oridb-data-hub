table oridbRefuted
"OriDB 3.0.0-trial.1 fragments that failed an ARS assay"
    (
    string  chrom;        "Reference sequence"
    uint    chromStart;   "Start, 0-based"
    uint    chromEnd;     "End, exclusive"
    string  name;         "Fragment identifier"
    uint    score;        "Not used"
    char[1] strand;       "Not applicable"
    uint    thickStart;   "Same as chromStart"
    uint    thickEnd;     "Same as chromEnd"
    uint    reserved;     "Colour, as itemRgb"
    string  result;       "Reported outcome of the assay"
    string  pmid;         "PubMed ID of the study"
    string  fragmentName; "Name used in the source study"
    string  otherNames;   "Other names for the same fragment"
    )
