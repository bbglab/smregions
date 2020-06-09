import numpy as np


def flatten_partitions(results):
    """Yield the partitions in the results"""
    for name, result in results.items():
        for partition in result['partitions']:
            yield (name, partition, result, np.random.randint(0, 2**32-1))


def partitions_list(total_size, chunk_size):
    """
    Create a list of values less or equal to chunk_size that sum total_size

    :param total_size: Total size
    :param chunk_size: Chunk size
    :return: list of integers
    """
    partitions = [chunk_size for _ in range(total_size // chunk_size)]

    res = total_size % chunk_size
    if res != 0:
        partitions += [res]

    return partitions


def compute_sampling(value):
    """Continue the computation from a partial chunck"""
    name, samples, result, seed = value

    muts_count = result['nmuts']
    items_to_simulate = result['simulation_items']
    items_to_simulate_prob = result['simulation_probs']
    regions_of_interest = result['region_of_interest']
    in_reg_counts = result['in_reg_counts']

    np.random.seed(seed)

    # Generate the simulated mutations
    indexes = range(len(items_to_simulate))
    background_index = np.random.choice(indexes, size=(samples, muts_count), p=items_to_simulate_prob,
                                        replace=True)

    # Get the mutations back from indexes
    list_mutations_simulated = []
    for list_muts in background_index:
        mutations = [items_to_simulate[x] for x in list_muts]
        list_mutations_simulated.append(mutations)

    for reg in regions_of_interest:
        name = reg.data
        start = reg.begin
        end = reg.end
        count_observed, count_simulated = in_reg_counts[name]
        for sim in list_mutations_simulated:
            count = 0
            for mut in sim:  # Get the mutations back from indexes
                if start <= mut[0] <= end:
                    count += 1
            count_simulated += count
        in_reg_counts[name] = (count_observed, count_simulated)

    return

