import math
import random
from otree.api import *


def compute_utility(player_choice, player_role, neighbors_choices):
    """
    compute player's utility as a function of
    - their choice (True=Blue, False=Red)
    - their role (minority=Blue, majority=Red)
    - their neighbors' choices (list of True/False)
    """
    from . import Constants

    num_neighbors = len(neighbors_choices)
    blue_neighbors = neighbors_choices.count(True)
    red_neighbors = neighbors_choices.count(False)

    #  minority player
    if player_role == "":
        return Constants.e if player_choice else 0

    # majority player
    if num_neighbors == 0:
        return Constants.s if not player_choice else 0

    p_blue = blue_neighbors / num_neighbors
    p_red = red_neighbors / num_neighbors

    if player_choice:  # Blue choice
        return Constants.z * (1 - math.exp(-Constants.lambda1 * p_blue)) / (1 - math.exp(-Constants.lambda1))
    else:  # Red choice
        return Constants.s + Constants.w * (1 - math.exp(-Constants.lambda2 * p_red)) / (1 - math.exp(-Constants.lambda2))

def payoff_table(degree):
    """
    Create a list of dictionaries showing z* and w* values
    for each possible number of coordinating neighbors.
    """
    from . import Constants

    if degree <= 0:
        return []

    table_data = []
    for n in range(degree + 1):
        p = n / degree
        zstar = Constants.z * (1 - math.exp(-Constants.lambda1 * p)) / (1 - math.exp(-Constants.lambda1))
        wstar = Constants.w * (1 - math.exp(-Constants.lambda2 * p)) / (1 - math.exp(-Constants.lambda2))
        table_data.append({
            'c_n': n,
            'zstar': round(zstar),
            'wstar': round(wstar),
        })
    return table_data
