version 1.0

struct Sublibrary {
    String sublibrary_id
	Array[File]+ fastq_R1s
	Array[File]+ fastq_R2s
}


struct ParseSample {
    String parse_sample_name
    File? vireo_donor_ids
}