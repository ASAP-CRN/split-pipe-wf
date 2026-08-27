version 1.0

# split-pipe workflow entrypoint

import "structs.wdl"
import "split_pipe/split_pipe.wdl" as Splitpipe

workflow split_pipe_wf {
  input {
    File input_R1
    File input_R2
  }

  call Splitpipe.split_pipe {
    input:
      input_R1 = input_R1,
      input_R2 = input_R2
  }

  output {
    File output_matrix = split_pipe.output_matrix
  }

  meta {
    description = "Main workflow for the split-pipe process"
  }

  parameter_meta {
    input_R1 = {
      description = "Input R1 FASTQ file"
    }
    input_R2 = {
      description = "Input R2 FASTQ file"
    }
  }

  runtime {
    docker: "splitpipe:1.8.2"
  }
}