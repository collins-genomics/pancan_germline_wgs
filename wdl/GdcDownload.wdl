# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Download a single file from NCI GDC using curl, then convert BAM to CRAM.
#
# 2026-09: GDC drops long-lived TLS connections (curl exit 56 mid-stream). The
# original `curl | samtools view -C` pipeline could not resume after a drop --
# a pipe has no on-disk byte offset to re-request from, and the downstream
# samtools cannot accept a re-spliced stream -- so every drop restarted
# 100+ GiB downloads from byte 0 and large files rarely finished.
# Now the BAM is landed on disk by a resuming curl loop (-C - re-requests from
# the current file size via an HTTP Range header), verified against the GDC
# index (size + md5), and only then converted. A download survives any number
# of connection drops as long as bytes keep advancing between attempts.


version 1.0

workflow GdcDownload {
  input {
    String uuid
    String filename
    File genomeRef
    File token
    String docker
    Int? diskGb
    Int preemptibleTries = 0
  }

  call DownloadBAM {
    input:
      uuid = uuid,
      filename = filename,
      token = token,
      genome_ref = genomeRef,
      docker = docker,
      diskGb = diskGb,
      preemptibleTries = preemptibleTries
  }

  output {
    File cram = DownloadBAM.cram
    File cram_idx = DownloadBAM.cram_idx
  }
}


task DownloadBAM {
  input {
    String uuid
    String filename
    File token
    File genome_ref
    String docker
    # BAM and CRAM coexist on disk during conversion; the largest BAM in the
    # TCGA cohorts is 274 GiB, so 400 was no longer safe.
    Int diskGb = 500
    # Downloads run for hours and a preempted VM loses the partial file
    # (resume only survives connection drops, not disk loss). Default to
    # on-demand; set >0 from the method config to use Spot VMs.
    Int preemptibleTries = 0
  }

  String cram_file = sub(filename, ".bam", ".cram")

  command <<<

    set -eu -o pipefail

    gdcTok=$(tr -d ' \t\r\n' < "~{token}")
    url='https://api.gdc.cancer.gov/data/~{uuid}'
    bam="~{filename}"

    # Expected size and md5 from the GDC index (public metadata, no token
    # required). Best effort: if unavailable, fall back to curl's own exit
    # status plus samtools quickcheck.
    meta=$(curl -sS --retry 5 --retry-delay 15 \
      "https://api.gdc.cancer.gov/files/~{uuid}?fields=file_size,md5sum" || true)
    want_size=$(printf '%s' "$meta" | grep -o '"file_size": *[0-9]*' | grep -o '[0-9]*' || true)
    want_md5=$(printf '%s' "$meta" | grep -o '"md5sum": *"[a-f0-9]\{32\}"' | grep -o '[a-f0-9]\{32\}' || true)
    echo "[gdc] uuid=~{uuid} expect size=${want_size:-?} md5=${want_md5:-?}" >&2

    # Resumable download. Only attempts that advance zero bytes count toward
    # failure; any forward progress resets the budget, so a flaky-but-moving
    # transfer is never abandoned.
    stall=0
    while :; do
      before=$(stat -c%s "$bam" 2>/dev/null || echo 0)
      if [ -n "$want_size" ] && [ "$before" -ge "$want_size" ]; then break; fi
      rc=0
      http=$(curl -sS -f -L -C - -o "$bam" -w '%{http_code}' \
              --connect-timeout 60 --speed-limit 102400 --speed-time 300 \
              -H "x-auth-token: $gdcTok" "$url") || rc=$?
      after=$(stat -c%s "$bam" 2>/dev/null || echo 0)
      echo "[gdc] $(date -u +%FT%TZ) curl rc=$rc http=${http:-?} bytes=$after/${want_size:-?}" >&2
      if [ "$rc" -eq 0 ] && { [ -z "$want_size" ] || [ "$after" -ge "$want_size" ]; }; then break; fi
      case "${http:-}" in
        401|403) echo "FATAL: GDC rejected the token (HTTP $http). Reissue it at portal.gdc.cancer.gov and replace the token file." >&2; exit 1 ;;
        404)     echo "FATAL: GDC has no file ~{uuid} (HTTP 404)" >&2; exit 1 ;;
      esac
      if [ "$rc" -eq 33 ] || [ "${http:-}" = 416 ]; then
        # Server refused the byte range: either the file is already complete,
        # or the partial is unusable and we must start over.
        if samtools quickcheck "$bam" 2>/dev/null; then break; fi
        echo "[gdc] resume refused and file incomplete; restarting from byte 0" >&2
        rm -f "$bam"
      fi
      if [ "$after" -gt "$before" ]; then stall=0; else stall=$((stall+1)); fi
      if [ "$stall" -ge 10 ]; then echo "FATAL: 10 consecutive attempts with no progress" >&2; exit 1; fi
      s=$((stall * 30)); [ "$s" -lt 5 ] && s=5; [ "$s" -gt 300 ] && s=300
      sleep "$s"
    done

    samtools quickcheck "$bam"
    if [ -n "$want_md5" ]; then
      got_md5=$(md5sum "$bam" | cut -c1-32)
      if [ "$got_md5" != "$want_md5" ]; then
        echo "FATAL: md5 mismatch after download: $got_md5 != $want_md5" >&2
        exit 1
      fi
    fi

    samtools view -@ 2 -x OQ --output-fmt-option level=6 -C -T "~{genome_ref}" -o "~{cram_file}" "$bam"

    rm -f "$bam"

    samtools index -@ 2 "~{cram_file}"

  >>>

  output {
    File cram = "~{cram_file}"
    File cram_idx = "~{cram_file + '.crai'}"
  }

  runtime {
    predefinedMachineType: "e2-standard-2"
    disks: "local-disk " + diskGb + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: preemptibleTries
    maxRetries: 3
  }
}
