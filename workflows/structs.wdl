version 1.0

struct Sublibrary {
    String sublibrary_id
	Array[File]+ fastq_R1s
	Array[File]+ fastq_R2s
}