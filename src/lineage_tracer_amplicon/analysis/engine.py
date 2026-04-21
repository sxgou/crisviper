"""
Main analysis engine for CRISPR lineage tracing amplicon analysis.
"""

import json
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path
import gzip

from ..core.structures import (
    AmpliconStructure, ReadResult, PrimerMatch, AnchorMatch,
    DeletionEvent, InsertionEvent, SNVEvent, Anchor, Target
)
from ..core.primer_validation import (
    validate_primers, create_primer_match
)
from ..core.anchor_alignment import (
    AnchorIndex, match_anchors, infer_deletions, find_closest_anchors
)
from ..core.mutation_classification import (
    classify_deletion_with_insertion, process_point_mutations,
    determine_target_state
)
from ..io.config_parser import load_config


class LineageTracerAnalyzer:
    """Main analyzer for CRISPR lineage tracing amplicon data."""
    
    def __init__(
        self,
        config_path: Path,
        primer_max_mismatch_rate: float = 0.15,
        min_anchor_identity: float = 0.8,
        kmer_size: int = 9,
        snv_window_size: int = 3,
        trim_primers: bool = False,  # Default False - do NOT trim
        min_deletion_size: int = 20,
        min_gap_size: int = 5      # Minimum gap to consider as deletion
    ):
        """
        Initialize analyzer with configuration.
        """
        self.structure = load_config(config_path)
        
        self.primer_max_mismatch_rate = primer_max_mismatch_rate
        self.min_anchor_identity = min_anchor_identity
        self.kmer_size = kmer_size
        self.snv_window_size = snv_window_size
        self.trim_primers = trim_primers
        self.min_deletion_size = min_deletion_size
        self.min_gap_size = min_gap_size
        
        self.primer5_seq, self.primer3_seq = self._extract_primers()
        self.anchors = self._extract_anchors()
        self.anchor_index = AnchorIndex(k=kmer_size)
        self.anchor_index.build(self.anchors)
        self.ref_length = len(self.structure.reference)
        
        self.stats = {
            'total_reads': 0,
            'passed_quality': 0,
            'primer5_detected': 0,
            'primer3_detected': 0,
            'both_primers_detected': 0,
            'passed_reads': 0,
            'snv_count': 0,
            'valid_snv_count': 0,
            'filtered_snv_count': 0,
            'deletion_count': 0,
            'insertion_count': 0,
            'per_target': {}
        }
    
    def _extract_primers(self) -> Tuple[str, str]:
        primer5_seq = ""
        primer3_seq = ""
        for feature in self.structure.features:
            if feature.get('type') == 'primer':
                name = feature.get('name', '')
                start = feature.get('start', 0)
                end = feature.get('end', 0)
                seq = self.structure.reference[start:end]
                if name == 'Primer5':
                    primer5_seq = seq
                elif name == 'Primer3':
                    primer3_seq = seq
        if not primer5_seq or not primer3_seq:
            raise ValueError("Primers not found in configuration")
        return primer5_seq, primer3_seq
    
    def _extract_anchors(self) -> List[Anchor]:
        anchors = []
        for feature in self.structure.features:
            feat_type = feature.get('type', '')
            name = feature.get('name', '')
            start = feature.get('start', 0)
            end = feature.get('end', 0)
            
            if feat_type in ['primer', 'prefix', 'postfix']:
                seq = self.structure.reference[start:end]
                anchor = Anchor(name=name, ref_start=start, ref_end=end, seq=seq, type=feat_type)
                anchors.append(anchor)
            elif feat_type == 'target':
                conserved = feature.get('conserved', [start, end])
                if len(conserved) == 2:
                    cons_start, cons_end = conserved
                    seq = self.structure.reference[cons_start:cons_end]
                    anchor = Anchor(name=f"{name}_conserved", ref_start=cons_start, ref_end=cons_end, seq=seq, type='conserved')
                    anchors.append(anchor)
            elif feat_type == 'pam_linker':
                pam = feature.get('pam', [start, end])
                if len(pam) == 2:
                    pam_start, pam_end = pam
                    seq = self.structure.reference[pam_start:pam_end]
                    anchor = Anchor(name=f"{name}_pam", ref_start=pam_start, ref_end=pam_end, seq=seq, type='pam')
                    anchors.append(anchor)
        return anchors
    
    def analyze_read(self, read_id: str, read_seq: str, quality: Optional[str] = None) -> ReadResult:
        self.stats['total_reads'] += 1
        
        # Primer validation
        has_p5, p5_match, has_p3, p3_match = validate_primers(
            read_seq, self.primer5_seq, self.primer3_seq,
            max_mismatch_rate=self.primer_max_mismatch_rate
        )
        
        if has_p5: self.stats['primer5_detected'] += 1
        if has_p3: self.stats['primer3_detected'] += 1
        
        is_valid = has_p5 and has_p3
        if is_valid:
            self.stats['both_primers_detected'] += 1
            self.stats['passed_reads'] += 1
        
        primer5_match_obj = create_primer_match('Primer5', p5_match, has_p5) if p5_match else None
        primer3_match_obj = create_primer_match('Primer3', p3_match, has_p3) if p3_match else None
        
        if not is_valid:
            return ReadResult(read_id=read_id, is_valid=False,
                             primer5_match=primer5_match_obj, primer3_match=primer3_match_obj)
        
        # Do NOT trim primers - use original read for anchor matching
        processed_seq = read_seq
        
        # Anchor matching
        anchor_matches = match_anchors(
            processed_seq, self.anchors, self.anchor_index,
            min_identity=self.min_anchor_identity, k=self.kmer_size
        )
        
        # Infer deletions with minimum gap size
        deletion_intervals = infer_deletions(
            anchor_matches, self.ref_length, min_gap_size=self.min_gap_size
        )
        
        # Classify deletion events
        deletion_events = []
        for del_start, del_end in deletion_intervals:
            left_anchor, right_anchor = find_closest_anchors(anchor_matches, del_start)
            read_fragment = ""
            if left_anchor and right_anchor:
                read_left = left_anchor.read_end
                read_right = right_anchor.read_start
                if read_right > read_left:
                    read_fragment = processed_seq[read_left:read_right]
            
            classification = classify_deletion_with_insertion(
                read_fragment, self.structure.reference, del_start, del_end
            )
            
            deletion_event = DeletionEvent(
                ref_start=del_start, ref_end=del_end,
                inserted_seq=classification['inserted_seq'],
                is_mmej=classification['is_mmej'],
                mh_left=classification['mh_left'], mh_right=classification['mh_right']
            )
            deletion_events.append(deletion_event)
            
            if classification['del_length'] >= self.min_deletion_size:
                self.stats['deletion_count'] += 1
        
        # Process SNVs
        anchor_match_info = []
        for match in anchor_matches:
            anchor_match_info.append({
                'read_start': match.read_start,
                'ref_start': match.anchor.ref_start,
                'length': match.read_end - match.read_start,
                'cigar': match.cigar
            })
        
        snv_events = process_point_mutations(
            processed_seq, self.structure.reference, anchor_match_info,
            self.structure, window_size=self.snv_window_size
        )
        self.stats['valid_snv_count'] += len(snv_events)
        
        # Determine target states
        target_states = {}
        for target in self.structure.targets:
            state = determine_target_state(
                target, deletion_events, [], snv_events  # insertion_events empty for now
            )
            target_states[target.name] = state
        
        return ReadResult(
            read_id=read_id, is_valid=True,
            primer5_match=primer5_match_obj, primer3_match=primer3_match_obj,
            anchor_matches=anchor_matches, deletion_events=deletion_events,
            insertion_events=[], snv_events=snv_events, target_states=target_states
        )
    
    def analyze_fastq(self, fastq_path: Path, output_dir: Path, max_reads: Optional[int] = None) -> Dict[str, Any]:
        output_dir.mkdir(parents=True, exist_ok=True)
        
        per_read_file = output_dir / "per_read_annotation.tsv"
        barcode_file = output_dir / "barcode_frequencies.tsv"
        stats_file = output_dir / "statistics.json"
        
        barcode_counts = {}
        edited_read_count = 0
        
        with open(per_read_file, 'w') as f_reads, open(barcode_file, 'w') as f_barcodes:
            # Headers
            f_reads.write("read_id\tvalid")
            for target in self.structure.targets:
                f_reads.write(f"\t{target.name}")
            f_reads.write("\tbarcode\tevents_json\n")
            f_barcodes.write("barcode\tcount\tfrequency\tsample\n")
            
            read_count = 0
            open_func = gzip.open if str(fastq_path).endswith('.gz') else open
            
            with open_func(fastq_path, 'rt') as f:
                for line1 in f:
                    if max_reads and read_count >= max_reads:
                        break
                    
                    read_id = line1.strip()[1:]
                    seq = f.readline().strip()
                    f.readline()  # +
                    quality = f.readline().strip()
                    
                    result = self.analyze_read(read_id, seq, quality)
                    
                    f_reads.write(f"{result.read_id}\t{result.is_valid}")
                    
                    if result.is_valid:
                        has_edit = False
                        for target in self.structure.targets:
                            state = result.target_states.get(target.name, 'WT')
                            if state != 'WT':
                                has_edit = True
                            f_reads.write(f"\t{state}")
                        
                        if has_edit:
                            edited_read_count += 1
                        
                        barcode_parts = [result.target_states.get(t.name, 'WT') for t in self.structure.targets]
                        barcode = '|'.join(barcode_parts)
                        barcode_counts[barcode] = barcode_counts.get(barcode, 0) + 1
                        
                        events_json = json.dumps({
                            'deletions': [{'start': d.ref_start, 'end': d.ref_end,
                                           'inserted_seq': d.inserted_seq, 'is_mmej': d.is_mmej}
                                          for d in result.deletion_events],
                            'snvs': [{'position': s.ref_pos, 'ref': s.ref_base,
                                      'alt': s.alt_base, 'target': s.target_name}
                                     for s in result.snv_events]
                        })
                        f_reads.write(f"\t{barcode}\t{events_json}\n")
                    else:
                        for _ in self.structure.targets:
                            f_reads.write("\t-")
                        f_reads.write("\t-\t-\n")
                    
                    read_count += 1
            
            total_valid = sum(barcode_counts.values())
            for barcode, count in barcode_counts.items():
                frequency = count / total_valid if total_valid > 0 else 0
                f_barcodes.write(f"{barcode}\t{count}\t{frequency:.6f}\tsample\n")
        
        if self.stats['passed_reads'] > 0:
            editing_efficiency = edited_read_count / self.stats['passed_reads']
        else:
            editing_efficiency = 0.0
        
        self.stats['editing_efficiency'] = editing_efficiency
        
        for target in self.structure.targets:
            self.stats['per_target'][target.name] = {'efficiency': 0.0, 'deletion': 0.0, 'insertion': 0.0, 'snv': 0.0}
        
        with open(stats_file, 'w') as f:
            json.dump(self.stats, f, indent=2)
        
        return self.stats