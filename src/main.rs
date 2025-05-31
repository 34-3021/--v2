use edlib_rs::edlibrs::*;
use std::fs::File;
use std::io::Read;

fn reverse_in_place(s: &mut String) {
    let mut chars: Vec<char> = s.chars().collect();
    chars.reverse();
    s.clear();
    s.extend(chars);
}

fn replace_dna_bases(input: &str) -> String {
    input
        .chars()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            _ => c,
        })
        .collect()
}

fn merge_intervals(intervals: &mut Vec<(usize, usize)>) -> (usize, usize) {
    if intervals.is_empty() {
        return (0, 0);
    }

    intervals.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    let mut merged = vec![];
    let mut largest_interval = (0, 0);

    let mut current = intervals[0];

    for &(start, end) in intervals.iter().skip(1) {
        if start <= current.1 {
            current.1 = current.1.max(end);
        } else {
            merged.push(current);
            if current.1 - current.0 > largest_interval.1 - largest_interval.0 {
                largest_interval = current;
            }
            current = (start, end);
        }
    }

    merged.push(current);
    if current.1 - current.0 > largest_interval.1 - largest_interval.0 {
        largest_interval = current;
    }

    largest_interval
}

fn paged_comparison(
    query_seq: &str,
    reference_seq_rc: &str,
    page_size: usize,
    threshold: i32,
) -> (
    Vec<(usize, usize)>,
    Vec<(usize, usize)>,
    (usize, usize),
    (usize, usize),
) {
    let len_ref = reference_seq_rc.len();
    let len_query = query_seq.len();

    let num_of_pages_ref = len_query.max(len_ref) / page_size + 1;
    let num_of_pages_query = len_query.max(len_query) / page_size + 1;

    let mut intervals_query = Vec::new();
    let mut intervals_ref = Vec::new();

    for i in 0..num_of_pages_query {
        let query = &query_seq[i * page_size..((i + 1) * page_size).min(len_query)];
        for j in 0..num_of_pages_ref {
            let target = &reference_seq_rc[j * page_size..((j + 1) * page_size).min(len_ref)];

            let align_res = edlibAlignRs(
                query.as_bytes(),
                target.as_bytes(),
                &EdlibAlignConfigRs::default(),
            );

            if align_res.editDistance < threshold {
                intervals_query.push((i * page_size, ((i + 1) * page_size).min(len_query)));
                intervals_ref.push((j * page_size, ((j + 1) * page_size).min(len_ref)));
            }
        }
    }

    let merged_interval_query = merge_intervals(&mut intervals_query);
    let merged_interval_ref = merge_intervals(&mut intervals_ref);

    (
        intervals_query,
        intervals_ref,
        merged_interval_query,
        merged_interval_ref,
    )
}

fn main() {
    let mut refer_file = File::open("data/reference.txt").unwrap();
    let mut refer_sequence = String::new();
    refer_file.read_to_string(&mut refer_sequence).unwrap();
    reverse_in_place(&mut refer_sequence);
    let reference_seq_rc = replace_dna_bases(&refer_sequence);

    let mut query_file = File::open("data/query.txt").unwrap();
    let mut query_seq = String::new();
    query_file.read_to_string(&mut query_seq).unwrap();

    let (_, __, merged_interval_query_1, merged_interval_ref_1) =
        paged_comparison(&query_seq, &reference_seq_rc, 128, 10);

    let query_result = merged_interval_query_1;
    let ref_result = (
        reference_seq_rc.len() - merged_interval_ref_1.1,
        reference_seq_rc.len() - merged_interval_ref_1.0
    );
    println!("({}, {}, {}, {})", query_result.0, query_result.1, ref_result.0, ref_result.1);
}