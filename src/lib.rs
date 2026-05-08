mod barcodes;
mod tags;

use pyo3::prelude::*;
use pyo3_polars::derive::polars_expr;
use pyo3_polars::export::polars_core::prelude::*;
use serde::Deserialize;

use tags::Tag;

// ---------------------------------------------------------------------------
// Kwargs schema
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct TagSpec {
    seq: String,
    info: String,
    max_mism: usize,
}

#[derive(Deserialize)]
struct FindTagsKwargs {
    tags: Vec<TagSpec>,
}

// ---------------------------------------------------------------------------
// Output type declaration
// ---------------------------------------------------------------------------

fn find_tags_output(_input_fields: &[Field]) -> PolarsResult<Field> {
    Ok(Field::new(
        "tags".into(),
        DataType::Struct(vec![
            Field::new("tag_seq".into(), DataType::String),
            Field::new("tag_type".into(), DataType::String),
            Field::new("tag_pos".into(), DataType::List(Box::new(DataType::UInt32))),
        ]),
    ))
}

// ---------------------------------------------------------------------------
// Expression
// ---------------------------------------------------------------------------

/// Find all tags in each sequence string and return a struct with three fields:
/// - tag_seq:  colon-joined matched sequences (e.g. "ACGT:TGCA")
/// - tag_type: colon-joined tag info strings  (e.g. "Dpm_r1:Term_r1")
/// - tag_pos:  list of start positions        (e.g. [0, 10])
#[polars_expr(output_type_func = find_tags_output)]
fn find_tags(inputs: &[Series], kwargs: FindTagsKwargs) -> PolarsResult<Series> {
    let tags: Vec<Tag> = kwargs
        .tags
        .iter()
        .map(|t| Tag::new(&t.seq, &t.info, t.max_mism))
        .collect();

    let ca = inputs[0].str()?;
    let len = ca.len();

    let mut tag_seq_builder = StringChunkedBuilder::new("tag_seq".into(), len);
    let mut tag_type_builder = StringChunkedBuilder::new("tag_type".into(), len);
    let mut pos_builder = ListPrimitiveChunkedBuilder::<UInt32Type>::new(
        "tag_pos".into(),
        len,
        len * 4,
        DataType::UInt32,
    );

    for opt in ca.iter() {
        match opt {
            None => {
                tag_seq_builder.append_value("");
                tag_type_builder.append_value("");
                pos_builder.append_slice(&[]);
            }
            Some(seq) => {
                let hits = barcodes::find_tags(seq.as_bytes(), &tags);
                if hits.is_empty() {
                    tag_seq_builder.append_value("");
                    tag_type_builder.append_value("");
                    pos_builder.append_slice(&[]);
                } else {
                    let tag_seq = hits
                        .iter()
                        .map(|(_, t)| t.get_seq())
                        .collect::<Vec<_>>()
                        .join(":");
                    let tag_type = hits
                        .iter()
                        .map(|(_, t)| t.get_info())
                        .collect::<Vec<_>>()
                        .join(":");
                    let positions: Vec<u32> = hits.iter().map(|(p, _)| *p as u32).collect();
                    tag_seq_builder.append_value(&tag_seq);
                    tag_type_builder.append_value(&tag_type);
                    pos_builder.append_slice(&positions);
                }
            }
        }
    }

    let tag_seq_s = tag_seq_builder.finish().into_series();
    let tag_type_s = tag_type_builder.finish().into_series();
    let tag_pos_s = pos_builder.finish().into_series();

    let out = StructChunked::from_series(
        "tags".into(),
        len,
        [&tag_seq_s, &tag_type_s, &tag_pos_s].iter().copied(),
    )?;

    Ok(out.into_series())
}

// ---------------------------------------------------------------------------
// Module registration
// ---------------------------------------------------------------------------

#[pymodule]
fn _lib(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
