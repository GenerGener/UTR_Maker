from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Dict, List, Tuple
import re

class UTRMaker:
    def __init__(self, genbank_file: str):
        """Initialize UTRMaker with a GenBank file."""
        self.record = next(SeqIO.parse(genbank_file, "genbank"))
        self.sequence = str(self.record.seq)
        self.features = self.record.features
        
    def find_segment_boundaries(self) -> Dict[str, List[Tuple[int, int]]]:
        """Find all locus segment boundaries in the GenBank file."""
        segments = {}
        
        for feature in self.features:
            if feature.type == "misc_feature":
                if "note" in feature.qualifiers:
                    note = feature.qualifiers["note"][0]
                    if "locus segment" in note.lower():
                        segment_num = re.search(r"locus segment (\d+)", note.lower())
                        if segment_num:
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            segment_id = f"segment_{segment_num.group(1)}"
                            
                            if segment_id not in segments:
                                segments[segment_id] = []
                            segments[segment_id].append((start, end))
                            
        return segments
    
    def find_coding_boundaries(self) -> Tuple[int, int]:
        """Find the first start and last stop of coding sequences."""
        starts = []
        stops = []
        
        for feature in self.features:
            if feature.type == "CDS":
                # Handle both simple and compound locations
                if hasattr(feature.location, 'parts'):
                    starts.append(int(feature.location.parts[0].start))
                    stops.append(int(feature.location.parts[-1].end))
                else:
                    starts.append(int(feature.location.start))
                    stops.append(int(feature.location.end))
        
        if starts and stops:
            return min(starts), max(stops)
        return None, None
    
    def extract_5prime_utr(self) -> SeqRecord:
        """Extract the 5' UTR sequence."""
        first_coding_start, _ = self.find_coding_boundaries()
        if first_coding_start is None:
            return None
            
        utr_seq = self.sequence[:first_coding_start]
        
        # Create a SeqRecord for the 5' UTR
        record = SeqRecord(
            Seq(utr_seq),
            id=f"{self.record.id}_5UTR",
            description=f"5' UTR from {self.record.id}"
        )
        return record
    
    def extract_3prime_utr(self) -> SeqRecord:
        """Extract the 3' UTR sequence."""
        _, last_coding_stop = self.find_coding_boundaries()
        if last_coding_stop is None:
            return None
            
        utr_seq = self.sequence[last_coding_stop:]
        
        # Create a SeqRecord for the 3' UTR
        record = SeqRecord(
            Seq(utr_seq),
            id=f"{self.record.id}_3UTR",
            description=f"3' UTR from {self.record.id}"
        )
        return record
    
    def save_utrs(self, output_prefix: str):
        """Save both UTRs to separate files."""
        # Extract and save 5' UTR
        utr5 = self.extract_5prime_utr()
        if utr5:
            with open(f"{output_prefix}_5UTR.fasta", "w") as f:
                SeqIO.write(utr5, f, "fasta")
                
        # Extract and save 3' UTR
        utr3 = self.extract_3prime_utr()
        if utr3:
            with open(f"{output_prefix}_3UTR.fasta", "w") as f:
                SeqIO.write(utr3, f, "fasta")
                
    def get_utr_details(self) -> Dict:
        """Get detailed information about the UTRs."""
        utr5 = self.extract_5prime_utr()
        utr3 = self.extract_3prime_utr()
        segments = self.find_segment_boundaries()
        
        return {
            "5_utr_length": len(utr5.seq) if utr5 else 0,
            "3_utr_length": len(utr3.seq) if utr3 else 0,
            "segments": segments
        }

# Example usage
if __name__ == "__main__":
    # Example with MZ242719
    maker1 = UTRMaker("MZ242719.txt")
    maker1.save_utrs("MZ242719")
    details1 = maker1.get_utr_details()
    print(f"MZ242719 5' UTR length: {details1['5_utr_length']}")
    print(f"MZ242719 3' UTR length: {details1['3_utr_length']}")
    
    # Example with MZ242720
    maker2 = UTRMaker("MZ242720.txt")
    maker2.save_utrs("MZ242720")
    details2 = maker2.get_utr_details()
    print(f"MZ242720 5' UTR length: {details2['5_utr_length']}")
    print(f"MZ242720 3' UTR length: {details2['3_utr_length']}")
