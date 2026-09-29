use anyhow::{Error, bail};
use rust_htslib::bam::{
    FetchDefinition, IndexedReader, Read, Writer as BamWriter,
    header::{Header, HeaderRecord},
    record::Record,
};

fn create_header() -> Header {
    let mut ret = Header::new();
    ret.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.4"));
    ret.push_comment(b"Generated with bamtax");
    ret
}

fn push_ref_to_header<T: ToString>(header: &mut Header, name: T, len: u64) {
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", name).push_tag(b"LN", len));
}

pub fn create_sub_bam<T: ToString + std::fmt::Display>(
    source: &mut IndexedReader,
    names: &[T],
    tids: &[u32],
    outname: &str,
) -> Result<(), Error> {
    let mut header = create_header();
    assert_eq!(names.len(), tids.len());

    for (name, tid) in std::iter::zip(names, tids) {
        let len = source.header().target_len(*tid).unwrap();
        push_ref_to_header(&mut header, name, len);
    }

    let mut out = BamWriter::from_path(outname, &header, rust_htslib::bam::Format::Sam)?;
    let mut rec = Record::new();
    for (i, tid) in tids.iter().enumerate() {
        source.fetch(FetchDefinition::Region(*tid as i32, 0, i64::MAX))?;

        loop {
            match source.read(&mut rec) {
                None => break,
                Some(Err(e)) => bail!("Failed to read from BAM: {e}"),
                Some(Ok(_)) => (),
            }

            rec.set_tid(i as i32);

            out.write(&rec)?;
        }
    }

    Ok(())
}
