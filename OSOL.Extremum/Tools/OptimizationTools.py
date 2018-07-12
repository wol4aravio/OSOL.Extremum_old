from Optimization.Algorithms.RandomSearch import RandomSearch
from Optimization.Algorithms.IntervalExplosionSearch import IntervalExplosionSearch
from Optimization.Algorithms.AdaptiveRandomSearch import AdaptiveRandomSearch
from Optimization.Algorithms.SimulatedAnnealing import SimulatedAnnealing
from Optimization.Algorithms.RandomSearchWithStatisticalAntiGradient import RandomSearchWithStatisticalAntiGradient
from Optimization.Algorithms.LuusJaakolaOptimization import LuusJaakolaOptimization
from Optimization.Algorithms.ModifiedHybridRandomSearch import ModifiedHybridRandomSearch


def create_algorithm_from_json(json_data):
    if 'RandomSearch' in json_data:
        return RandomSearch.from_json(json_data)
    elif 'IntervalExplosionSearch' in json_data:
        return IntervalExplosionSearch.from_json(json_data)
    elif 'AdaptiveRandomSearch' in json_data:
        return AdaptiveRandomSearch.from_json(json_data)
    elif 'SimulatedAnnealing' in json_data:
        return SimulatedAnnealing.from_json(json_data)
    elif 'RandomSearchWithStatisticalAntiGradient' in json_data:
        return RandomSearchWithStatisticalAntiGradient.from_json(json_data)
    elif 'LuusJaakolaOptimization' in json_data:
        return LuusJaakolaOptimization.from_json(json_data)
    elif 'ModifiedHybridRandomSearch' in json_data:
        return ModifiedHybridRandomSearch.from_json(json_data)
    else:
        raise Exception('Unsupported Optimization Algorithm')
