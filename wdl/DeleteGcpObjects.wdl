# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Parallel deletion of objects stored in GCP buckets


version 1.0


workflow DeleteGcpObjects {
  input {
    File uri_list
    Int n_cpu = 4
    Int uris_per_shard = 1000000
    String linux_docker = "ubuntu:latest"
    String sdk_docker = "google/cloud-sdk:latest"
  }

  call ShuffleAndShard {
    input:
      uri_list = uri_list,
      uris_per_shard = uris_per_shard,
      linux_docker = linux_docker
  }

  scatter ( shard in ShuffleAndShard.uri_list_shards ) {
    call DeleteObjects {
      input:
        uri_list = shard,
        n_cpu = n_cpu,
        sdk_docker = sdk_docker
    }
  }
}


task ShuffleAndShard {
  input {
    File uri_list
    Int uris_per_shard
    String linux_docker
  }

  Int input_gb = ceil(size(uri_list, "GB"))
  Int disk_gb = 10 + ceil(3 * input_gb)
  Float mem_gb = 2 + ceil(3 * input_gb)
  Int n_cpu = ceil(mem_gb / 4)

  command <<<
    set -eu -o pipefail

    shuf ~{uri_list} \
    | split -d -a 8 -l ~{uris_per_shard} - sharded_uris_
  >>>

  output {
    Array[File] uri_list_shards = glob("sharded_uris_*")
  }

  runtime {
    docker: linux_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 5
  }
}


task DeleteObjects {
  input {
    File uri_list
    Int n_cpu
    String sdk_docker
  }

  Int n_threads = 2 * n_cpu
  Float mem_gb = 2 + ( 0.5 * n_cpu )
  Int disk_gb = 10 + ceil(2 * size(uri_list, "GB"))

  command <<<
    cat ~{uri_list} \
    | gsutil \
        -m \
        -o GSUtil:parallel_process_count=1 \
        -o GSUtil:parallel_thread_count=~{n_threads} \
        rm -fI
  >>>

  runtime {
    docker: sdk_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 5
  }
}

