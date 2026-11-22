table oridbElements
"OriDB 3.0.0-trial.1 origin elements"
    (
    string  chrom;        "Reference sequence"
    uint    chromStart;   "Start, 0-based"
    uint    chromEnd;     "End, exclusive"
    string  name;         "OriDB element accession"
    uint    score;        "Motif score as a percentage of the matrix maximum, times 10"
    char[1] strand;       "Strand the element matches on"
    uint    thickStart;   "Same as chromStart"
    uint    thickEnd;     "Same as chromEnd"
    uint    reserved;     "Element-class colour, as itemRgb"
    string  elementClass; "A, B2 or B3"
    string  evidence;     "confirmed or predicted"
    string  origin;       "Accession of the origin this element belongs to"
    string  studyPmid;    "PubMed ID of the study that defined it"
    string  motifScore;   "Motif score as a percentage, or '.' where not scored"
    )
