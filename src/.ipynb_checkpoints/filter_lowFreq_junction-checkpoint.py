import pandas as pd
import multiprocessing as mp
from typing import List, Set
from collections import defaultdict


class FilterLowfreqJunction:
    """
    Filters junctions with frequency below a specified threshold.
    
    Args:
        min_junction_freq: Minimum frequency threshold for keeping junctions
        num_processes: Number of processes for parallel processing
    """
    
    def __init__(self, min_junction_freq: int = 2, num_processes: int = 10):
        self.min_junction_freq = min_junction_freq
        self.num_processes = num_processes
    
    @staticmethod
    def exon_to_junction(exon_str: str) -> List[str]:
        """Convert exon chain string to list of junctions."""
        exons = exon_str.split('-')
        return [f"{exons[i]}-{exons[i+1]}" for i in range(0, len(exons)-1, 2)]
    
    def _get_low_freq_junctions(self, df_chr: pd.DataFrame) -> Set[str]:
        """Identify low-frequency junctions for a chromosome."""
        exploded = df_chr.explode('junctionChain')
        junction_counts = exploded.groupby('junctionChain')['frequency'].sum()
        return set(junction_counts[junction_counts <= self.min_junction_freq].index)
    
    def _filter_chromosome(self, df_chr: pd.DataFrame) -> pd.DataFrame:
        """Filter low-frequency junctions for a single chromosome."""
        low_freq_junctions = self._get_low_freq_junctions(df_chr)
        
        # Create mask for rows to keep (inverse of rows with low-freq junctions)
        mask = ~df_chr['junctionChain'].apply(
            lambda junctions: any(j in low_freq_junctions for j in junctions)
        )
        
        return df_chr[mask]
    
    def filterJun(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Main filtering method that processes the entire DataFrame.
        
        Args:
            df: Input DataFrame containing exonChain and frequency columns
            
        Returns:
            Filtered DataFrame with low-frequency junctions removed
        """
        # Convert exon chains to junctions if not already done
        if 'junctionChain' not in df.columns:
            df = df.assign(junctionChain=df['exonChain'].apply(self.exon_to_junction))
        
        # Split by chromosome for parallel processing
        chromosomes = df['Chr'].unique()
        grouped = [df[df['Chr'] == chr_] for chr_ in chromosomes]
        
        # Process chromosomes in parallel
        with mp.Pool(self.num_processes) as pool:
            results = pool.map(self._filter_chromosome, grouped)
            
        return pd.concat(results, ignore_index=True)