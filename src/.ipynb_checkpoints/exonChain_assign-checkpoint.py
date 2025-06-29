from typing import List, Tuple, Dict
from collections import defaultdict
import itertools

class ExonChainAssign:
    def __init__(self, delta: int = 10):
        self.delta = delta
        self.index: Dict[Tuple[int, ...], List[str]] = defaultdict(list)
        self.chain_set = set()

    def parse_chain(self, chain: str) -> List[int]:
        return list(map(int, chain.split("-")))

    def get_intervals(self, positions: List[int]) -> List[Tuple[int, int]]:
        return [(positions[i], positions[i+1]) for i in range(len(positions) - 1)]

    def to_bucket(self, s: int, e: int) -> Tuple[int, int]:
        if self.delta == 0:
            return (s, e)  # 完全匹配时直接用原始值，不做bucket划分
        return (s // self.delta, e // self.delta)

    def get_bucketed_subchains(self, chain: str) -> List[Tuple[Tuple[int, int], ...]]:
        intervals = self.get_intervals(self.parse_chain(chain))
        buckets = [self.to_bucket(s, e) for s, e in intervals]
        return [tuple(buckets[i:i+1]) for i in range(len(buckets) - 1 + 1)]

    def build_index(self, chains: List[str]):
        for chain in chains:
            for key in self.get_bucketed_subchains(chain):
                self.index[key].append(chain)
                
    def build_index_fortrun(self, chains: list):
        self.chain_set = set(chains)

    def get_fuzzy_keys(self, key: Tuple[Tuple[int, int], ...]) -> List[Tuple[Tuple[int, int], ...]]:
        options = []
        for s, e in key:
            options.append([
                (s - 1, e - 1), (s - 1, e), (s - 1, e + 1),
                (s, e - 1),     (s, e),     (s, e + 1),
                (s + 1, e - 1), (s + 1, e), (s + 1, e + 1)
            ])
        return list(itertools.product(*options))

    def is_fuzzy_assign(self, query: List[Tuple[int, int]], target: List[Tuple[int, int]]) -> bool:
        if len(query) != len(target):
            return False
        for (s1, e1), (s2, e2) in zip(query, target):
            if abs(s1 - s2) > self.delta or abs(e1 - e2) > self.delta:
                return False
        return True

    def assign_exonChain(self, query_chain: str) -> List[str]:
        query_intervals = self.get_intervals(self.parse_chain(query_chain))
        query_buckets = [self.to_bucket(s, e) for s, e in query_intervals]

        candidates = set()
        if len(query_buckets) < 1:
            return []

        for i in range(len(query_buckets) - 1 + 1):
            sub = tuple(query_buckets[i:i+1])
            for fuzzy_key in self.get_fuzzy_keys(sub):
                for cand in self.index.get(fuzzy_key, []):
                    candidates.add(cand)

        assignes = []
        for cand in candidates:
            cand_intervals = self.get_intervals(self.parse_chain(cand))
            if len(cand_intervals) != len(query_intervals):
                continue
            if self.is_fuzzy_assign(query_intervals, cand_intervals):
                assignes.append(cand)

        return assignes
    
    def assign_truncation(self, full_chain: str) -> list:
        positions = self.parse_chain(full_chain)
        assignes = []

        # 提取所有长度 >= 2 的连续子链
        for i in range(len(positions) - 1):
            for j in range(i + 2, len(positions) + 1):  # j 是 end+1
                subchain = "-".join(map(str, positions[i:j]))
                if subchain in self.chain_set:
                    assignes.append(subchain)

        return assignes