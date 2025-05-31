use core::panic;
use std::collections::HashMap;
use std::io;

fn reverse_dna(dna: &str) -> String {
    dna.chars().rev().map(|c| match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => panic!("Invalid DNA character")
    }).collect()
}

fn dna_to_num(dna: char) -> u64 {
    match dna {
        'A' => 1,
        'T' => 2,
        'C' => 3,
        'G' => 4,
        _ => panic!("Invalid DNA character")
    }
}

#[derive(Debug, Clone)]
struct RefSeq {
    start: u64,
    end: u64,
    reverse: bool
}

const MOD: u64 = 1_000_000_000_000_7;

fn insert_subseq(dna: &str, map: &mut HashMap<u64, RefSeq>, reverse: bool) {
    let dna_len = dna.len();
    let dna: Vec<char> = if reverse {
        reverse_dna(dna)
    } else {
        dna.to_string()
    }.chars().collect();

    for start in 0..dna_len {
        let mut hash: u64 = 0;
        for end in start..dna_len {
            hash = (hash * 5 + dna_to_num(dna[end])) % MOD;
            map.entry(hash).or_insert(RefSeq {
                start: if reverse { dna_len - end - 1 } else { start } as u64,
                end: if reverse { dna_len - start - 1 } else { end } as u64,
                reverse
            });
        }
    }
}

#[derive(Debug, Clone)]
struct Trace {
    ref_seq: RefSeq,
    next: u64,
}

fn search_trace(query: &str, ref_map: &HashMap<u64, RefSeq>) -> Vec<Option<Trace>> {
    let query_len = query.len();
    let query: Vec<char> = query.chars().collect();
    let mut dp = vec![u64::MAX - 20; query_len + 1]; // -20 to avoid +1 overflow
    dp[query_len] = 0;

    let mut trace = vec![None; query_len + 1];
    for start in (0..query_len).rev() {
        let mut hash: u64 = 0;
        for end in start..query_len {
            hash = (hash * 5 + dna_to_num(query[end])) % MOD;
            if let Some(ref_seq) = ref_map.get(&hash) {
                if (dp[start] > dp[end + 1] + 1) || (dp[start] == dp[end + 1] + 1 && !ref_seq.reverse) {
                    dp[start] = dp[end + 1] + 1;
                    trace[start] = Some(Trace {
                        ref_seq: ref_seq.clone(),
                        next: (end + 1) as u64
                    });
                }
            }
        }
    }

    trace
}

fn parse_trace(trace: Vec<Option<Trace>>, query_len: usize) -> Vec<RefSeq> {
    let mut result = Vec::new();
    let mut pos = 0;
    while pos < query_len {
        match &trace[pos] {
            None => panic!("No correct trace found"),
            Some(Trace { ref_seq, next }) => {
                result.push(ref_seq.clone());
                pos = *next as usize;
            }
        }
    }
    result
}

fn main() {
    let mut lines = io::stdin().lines().map(|l| l.unwrap().trim().to_string());
    let ref_ = lines.next().unwrap();
    let query = lines.next().unwrap();

    let mut ref_map = HashMap::new();
    insert_subseq(&ref_, &mut ref_map, false);
    insert_subseq(&ref_, &mut ref_map, true);

    let trace = search_trace(&query, &ref_map);
    let result = parse_trace(trace, query.len());

    result.iter().for_each(|x| {
        if x.start == x.end {
            println!("Matched sequence: ref[{}]=\"{}\", reverse={}", x.start, ref_.chars().nth(x.start as usize).unwrap(), x.reverse);
        } else {
            println!("Matched sequence: ref[{}~{}]=\"{}\", reverse={}", x.start, x.end, ref_.get(x.start as usize ..= x.end as usize).unwrap(), x.reverse);
        }
    });
}
