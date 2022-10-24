import logging

import numpy as np

from smregions import reference
from smregions.walker import partitions_list
from collections import defaultdict

logger = logging.getLogger(__name__)


class ElementExecutor:
    """
    Executors that do the analysis per genomic element.

    Args:
        element_id (str): element ID
        muts (list): list of mutations belonging to that element (see :ref:`mutations <mutations dict>` inner items)
        segments (list): list of segments belonging to that element (see :ref:`elements <elements dict>` inner items)
        signature (dict): probabilities of each mutation (see :ref:`signatures <signature dict>`)
        config (dict): configuration

    """

    def __init__(self, element_id, muts, segments, regions_of_interest, signature, config, mutpersample, seed=None):
        # Input attributes
        self.name = element_id

        self.muts = muts
        self.signature = signature
        self.segments = segments
        self.regions_of_interest = regions_of_interest

        # Configuration parameters
        self.mutpersample = mutpersample
        self.sampling_size = config['sampling']
        self.sampling_chunk = config['sampling_chunk'] * 10**6
        self.seed = seed

        # Output attributes
        self.result = {}

    def run(self):

        np.random.seed(self.seed)

        observed = [(m['POSITION'], m['ALT']) for m in self.muts]

        muts_count = len(self.muts)
        self.result['nmuts'] = muts_count
        self.result['partitions'] = []
        self.result['in_reg_counts'] = {}
        if muts_count > 0 and len(self.regions_of_interest) > 0:

            items_to_simulate = []
            items_to_simulate_prob = defaultdict(list)
            in_reg_counts = {}

            for segment in self.segments:
                for pos, ref_triplet in zip(range(segment['START'], segment['END'] + 1),
                                            reference.generate_triplets(segment['CHROMOSOME'], segment['START'], segment['END'])):

                    if 'N' in ref_triplet:
                        continue
                    ref = ref_triplet[1]

                    for alt in'ACGT':
                        if alt == ref:
                            continue

                        items_to_simulate.append((pos, alt))

                        #mEdit +. Instead a list of simulations-probabilities, now we have a list of list (same but per sample)
                        for sample_id in self.mutpersample.keys():
                            if self.signature[sample_id] is None:
                                items_to_simulate_prob[sample_id].append(1.0)
                            else:
                                items_to_simulate_prob[sample_id].append(self.signature[sample_id].get(ref_triplet + '>' + alt, 0.0))

            for sample_id in self.mutpersample.keys():
                if sum(items_to_simulate_prob[sample_id]) == 0:
                    logger.warning('Probability of simulation equal to 0 in {}'.format(self.name))
                else:
                    list_mutations_simulated = []
                    for sample_id, sample_mutations in self.mutpersample.items():
                        # normalize probs
                        items_to_simulate_prob[sample_id] = np.array(items_to_simulate_prob[sample_id])
                        items_to_simulate_prob[sample_id] = items_to_simulate_prob[sample_id] / sum(items_to_simulate_prob[sample_id])

                        # Calculate sampling parallelization partitions
                        chunk_count = (self.sampling_size * muts_count) // self.sampling_chunk
                        chunk_size = self.sampling_size if chunk_count == 0 else self.sampling_size // chunk_count
                        self.result['partitions'] = partitions_list(self.sampling_size, chunk_size)

                        # Run first partition
                        first_partition = self.result['partitions'].pop(0)

                        indexes = range(len(items_to_simulate))
                        background_index = np.random.choice(indexes, size=(first_partition, len(sample_mutations)), p=items_to_simulate_prob[sample_id], replace=True)

                        # Iterate over the regions of that particular element

                        # Get the mutations back from indexes
                        for list_muts in background_index:
                            mutations = [items_to_simulate[x] for x in list_muts]
                            list_mutations_simulated.append(mutations)

                    for reg in self.regions_of_interest:
                        name = reg.data
                        start = reg.begin
                        end = reg.end
                        count_observed, count_simulated = 0, 0
                        for mut in observed:
                            if start <= mut[0] <= end:
                                count_observed += 1
                        for sim in list_mutations_simulated:
                            count = 0
                            for mut in sim:  # Get the mutations back from indexes
                                if start <= mut[0] <= end:
                                    count += 1
                            count_simulated += count
                        in_reg_counts[name] = (count_observed, count_simulated)

                    # Sampling parallelization (if more than one partition)
                    if len(self.result['partitions']) > 0:
                        self.result['region_of_interest'] = self.regions_of_interest
                        self.result['simulation_items'] = items_to_simulate
                        self.result['simulation_probs'] = items_to_simulate_prob

                self.result['in_reg_counts'] = in_reg_counts

        return self
