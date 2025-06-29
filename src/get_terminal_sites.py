import pandas as pd
import numpy as np
import multiprocessing as mp
from sklearn.cluster import DBSCAN
from functools import partial
from scipy.spatial.distance import cdist
from collections import Counter

class TerminalSitesProcessor:
    def __init__(self, cluster_group_size=1500, eps=10, min_samples=20, num_processes=10):
        self.cluster_group_size = cluster_group_size
        self.eps = eps
        self.min_samples = min_samples
        self.num_processes = num_processes

    def get_representative_sites_DBSCAN(self, values, mode='min'):
        if len(values) >= self.cluster_group_size:
            values = values.sample(n=self.cluster_group_size, random_state=42)

        X = np.array(values).reshape(-1, 1)
        min_samples = max(int(len(X) * 0.1), self.min_samples)
        labels = DBSCAN(eps=self.eps, min_samples=min_samples).fit_predict(X)

        df = pd.DataFrame({'value': values, 'cluster': labels})
        valid_clusters = [c for c in df['cluster'].unique() if c != -1]

        if not valid_clusters:
            return [df['value'].mode().min()] if mode == 'min' else [df['value'].mode().max()]
        elif len(valid_clusters) == 1:
            modes = df[df['cluster'] != -1]['value'].mode()
            return [modes.mode().min()] if mode == 'min' else [modes.mode().max()]
        else:
            results = []
            for c in valid_clusters:
                mode_vals = df[df['cluster'] == c]['value'].mode()
                result = mode_vals.min() if mode == 'min' else mode_vals.max()
                results.append(result)
            return sorted(set(results))

    def get_terminal_sites_forChr(self, df):
        results = []

        for _, row in df.iterrows():
            starts = pd.Series(row['TrStart_reads'])
            ends = pd.Series(row['TrEnd_reads'])

            pred_starts = self.get_representative_sites_DBSCAN(starts, mode='min')
            pred_ends = self.get_representative_sites_DBSCAN(ends, mode='max')

            raw_pairs = set(zip(row['TrStart_reads'], row['TrEnd_reads']))

            if len(pred_starts) == 1 and len(pred_ends) == 1:
                s, e = pred_starts[0], pred_ends[0]
                results.append({
                    'Chr': row['Chr'],
                    'Strand': row['Strand'],
                    'SSC': row['SSC'],
                    'Group': row['Group'],
                    'TrStart': s,
                    'TrEnd': e,
                    'frequency': row['frequency']
                })
                continue

            valid_pairs = [(s, e) for s in pred_starts for e in pred_ends if (s, e) in raw_pairs]

            if len(valid_pairs) == 1:
                s, e = valid_pairs[0]
                results.append({
                    'Chr': row['Chr'],
                    'Strand': row['Strand'],
                    'SSC': row['SSC'],
                    'Group': row['Group'],
                    'TrStart': s,
                    'TrEnd': e,
                    'frequency': row['frequency']
                })
            elif len(valid_pairs) > 1:
                raw_points = np.array(list(raw_pairs))
                centers = np.array(valid_pairs)

                dist_matrix = cdist(raw_points, centers, metric='euclidean')
                assigned_centers = dist_matrix.argmin(axis=1)

                center_counter = Counter(assigned_centers)
                total = sum(center_counter.values())

                for idx, (s, e) in enumerate(valid_pairs):
                    freq = row['frequency'] * (center_counter[idx] / total) if idx in center_counter else 0
                    results.append({
                        'Chr': row['Chr'],
                        'Strand': row['Strand'],
                        'SSC': row['SSC'],
                        'Group': row['Group'],
                        'TrStart': s,
                        'TrEnd': e,
                        'frequency': freq
                    })
            else:
                s = min(pred_starts)
                e = max(pred_ends)
                results.append({
                    'Chr': row['Chr'],
                    'Strand': row['Strand'],
                    'SSC': row['SSC'],
                    'Group': row['Group'],
                    'TrStart': s,
                    'TrEnd': e,
                    'frequency': row['frequency']
                })

        return pd.DataFrame(results)

    def get_terminal_sites(self, df):
        dfChr_list = [dfChr for _, dfChr in df.groupby(['Chr','Strand'], observed=True)]
        partial_func = partial(self.get_terminal_sites_forChr)

        with mp.Pool(self.num_processes) as pool:
            results = pool.map(partial_func, dfChr_list)

        df_result = pd.concat(results, ignore_index=True).reset_index(drop=True)
        
        return df_result
