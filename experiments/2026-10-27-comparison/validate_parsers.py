"""Cross-check native-format readers against independent text output/read IDs."""
import csv
import itertools
import json
from pathlib import Path
import struct
import subprocess

from make_inputs import DATA, RAW, ROOT, tag, MASK
from read_outputs import constellation_records, pack_dna, rad_records, unpack_dna


def validate_constellation():
    states={1:"unique_gene",2:"ambiguous_gene",3:"ambiguous_transcript_same_gene",
            4:"antisense_gene",5:"unmapped",6:"low_complexity",7:"low_quality"}
    old=ROOT/"benchdata/throughput-2026-10-27"
    checked=0
    with (old/"optimized.tsv").open() as f:
        for row,record in itertools.islice(zip(csv.DictReader(f,delimiter="\t"),
                 constellation_records(old/"optimized.rad")),100_000):
            assert states[record["state"]]==row["assignment_type"]
            assert record["primary"]==(int(row["gene_id"]) if row["gene_id"]!="." else 0xffffffff)
            assert record["score"]==int(row["score"]) and record["flags"]==int(row["flags"])
            if record["valid"]==(1<<28)-1:
                assert unpack_dna(record["bc"],16)==row["cell_barcode"]
                assert unpack_dna(record["umi"],12)==row["umi"]
            checked+=1
    assert checked==100_000
    return checked


def tag_desc(name,kind):
    return struct.pack("<H",len(name))+name.encode()+bytes([kind])


def validate_piscem():
    source=DATA/"synthetic-piscem/map.rad"
    metadata={}
    iterator=rad_records(source,metadata)
    first=next(iterator)
    assert metadata["version"]==2
    # Copy the original reference header and first chunk WITHOUT altering record
    # bytes; remove only the new schema prelude/roles so the installed legacy
    # alevin-fry reader can independently decode the same payload.
    with source.open("rb") as f:
        f.seek(metadata["header_start"])
        reference_header=f.read(metadata["chunk_count_offset"]-metadata["header_start"])
        f.seek(metadata["chunks_offset"])
        chunk_header=f.read(8)
        chunk_size,nrecords=struct.unpack("<II",chunk_header)
        payload=f.read(chunk_size-8)
    destination=DATA/"validation-first-chunk-legacy.rad"
    with destination.open("wb") as f:
        f.write(reference_header+struct.pack("<QH",1,2))
        f.write(tag_desc("cblen",2)+tag_desc("ulen",2))
        f.write(struct.pack("<H",2)+tag_desc("b",3)+tag_desc("u",3))
        f.write(struct.pack("<H",1)+tag_desc("compressed_ori_refid",3))
        f.write(struct.pack("<HH",16,12)+chunk_header+payload)
    expected=[first,*itertools.islice(iterator,nrecords-1)]
    iterator.close()
    assert len(expected)==nrecords
    assert list(rad_records(destination))==expected
    cmd=[str(ROOT/".conda-af/bin/alevin-fry"),"view","--rad",str(destination)]
    hit_count=0
    with (RAW/"parser-alevin-view.log").open("w") as log:
        with subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log,text=True) as process:
            for line in process.stdout:
                if not line.startswith("ID:"):
                    continue
                parts=line.rstrip().split("\t")
                fields=dict(part.split(":",1) for part in parts[:-1])
                row=int(fields["ID"]); hit=int(fields["HI"])-1
                bc,umi,hits=expected[row]
                assert fields["CB"]==unpack_dna(bc,16)
                assert fields["UMI"]==unpack_dna(umi,12)
                assert int(fields["NH"])==len(hits)
                assert (parts[-1],fields["DIR"]=="true")==hits[hit]
                hit_count+=1
            assert process.wait()==0
    assert hit_count==sum(len(hits) for bc,umi,hits in expected)
    return dict(records=nrecords,alignments=hit_count,legacy_reader=str(cmd[0]))


if __name__=="__main__":
    assert all(pack_dna(tag(i))^MASK==i for i in range(1_000_000))
    report=dict(identity_roundtrips=1_000_000,
                constellation_records_checked_against_tsv=validate_constellation(),
                piscem_checked_against_alevin_view=validate_piscem(),
                kallisto_identity_check="Every BUS flag (-n) is checked against the decoded UMI in read_outputs.py")
    (RAW/"parser-validation.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report),flush=True)
