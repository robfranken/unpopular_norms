from otree.api import *

#custom functions:
import math
import random

def compute_utility(player_choice, player_role, neighbors_choices, Constants):

    num_neighbors = len(neighbors_choices)
    blue_neighbors = neighbors_choices.count(True)
    red_neighbors = neighbors_choices.count(False)

    # minority player
    if player_role == Constants.minority:
        return Constants.e if player_choice else 0

    # majority player
    if num_neighbors == 0:
        return Constants.s if not player_choice else 0

    p_blue = blue_neighbors / num_neighbors
    p_red = red_neighbors / num_neighbors

    if player_choice:
        return Constants.z * (1 - math.exp(-Constants.lambda1 * p_blue)) / (
            1 - math.exp(-Constants.lambda1)
        )
    else:
        return Constants.s + Constants.w * (1 - math.exp(-Constants.lambda2 * p_red)) / (
            1 - math.exp(-Constants.lambda2)
        )


def payoff_table(degree, Constants):

    if degree <= 0:
        return []

    table_data = []

    for n in range(degree + 1):

        p = n / degree

        zstar = Constants.z * (1 - math.exp(-Constants.lambda1 * p)) / (
            1 - math.exp(-Constants.lambda1)
        )

        wstar = Constants.w * (1 - math.exp(-Constants.lambda2 * p)) / (
            1 - math.exp(-Constants.lambda2)
        )

        table_data.append({
            "c_n": n,
            "zstar": round(zstar),
            "wstar": round(wstar),
        })

    return table_data



doc = """
They receive a brief (role-based) instruction, after which they complete a set of comprehension questions.
"""

class Constants(BaseConstants):
    title = "The Fashion Dilemma"
    name_in_url = 'comprehension'
    players_per_group = None
    num_rounds = 1
    majority = "Red"
    minority = "Blue"
    s = 15
    e = 10
    z = 50
    w = 40
    lambda1 = 4.3
    lambda2 = 1.8
    base_payment = 2.50
    max_payment = 7.50
    points_per_euro_majority = 167
    points_per_euro_minority = 33
    decision_pages_timeout_seconds = 20
    other_pages_timeout_seconds = 15
    introduction_timeout_seconds = 10*60
    comprehension_timeout_seconds = 5*60
    max_retries = 3
    nrounds = 25

class Subsession(BaseSubsession):
    pass

class Group(BaseGroup):
    pass

class Player(BasePlayer):
    q_red_zero = models.IntegerField(min=0, label="")
    q_blue_zero = models.IntegerField(min=0, label="")
    q_red_half = models.IntegerField(min=0, label="")
    q_blue_half = models.IntegerField(min=0, label="")

    payoff_red_zero = models.IntegerField()
    payoff_blue_zero = models.IntegerField()
    payoff_red_half = models.IntegerField()
    payoff_blue_half = models.IntegerField()

    #we also count the number of wrong submissions during the comprehension check (increments per wrong submission)
    comprehension_retries = models.IntegerField(initial=0)

class IntroductionPage(Page):

    #timeout_seconds = Constants.introduction_timeout_seconds

    @staticmethod
    def is_displayed(player):
        return player.participant.vars.get('consent')

    @staticmethod
    def vars_for_template(player):
        degree = 2 # for instruction, assume 2 neighbors (this can be tweaked)
        table_data = payoff_table(degree, Constants)

        return dict(
            role=player.participant.role,
            group_size=50,
            degree=degree,
            range_neighbors=list(range(degree + 1)) if degree > 0 else [],
            table_data=table_data,
            base="{:.2f}".format(Constants.base_payment),
            max="{:.2f}".format(Constants.max_payment),
            num_rounds_lower=round(Constants.nrounds * 0.9),
            num_rounds_upper = round(Constants.nrounds * 1.1),
        )

class ComprehensionPage(Page):
    form_model = 'player'
    form_fields = ['q_red_zero', 'q_blue_zero', 'q_red_half', 'q_blue_half']

    @staticmethod
    def is_displayed(player):
        return player.participant.vars.get('consent')

    #timeout_seconds = Constants.comprehension_timeout_seconds

    #@staticmethod
    #def get_timeout_seconds(player):
    #    if player.comprehension_retries >= Constants.max_retries:
    #        return 1
    #    return Constants.comprehension_timeout_seconds

    def vars_for_template(player):
        degree = 2
        table_data = payoff_table(degree, Constants)

        neighbors_all_blue = [True] * degree
        neighbors_half_half = [True] * (degree // 2) + [False] * (degree - degree // 2)

        role = player.participant.role

        payoff_red_zero = compute_utility(False, role, neighbors_all_blue, Constants)
        payoff_blue_zero = compute_utility(True, role, neighbors_all_blue, Constants)
        payoff_red_half = compute_utility(False, role, neighbors_half_half, Constants)
        payoff_blue_half = compute_utility(True, role, neighbors_half_half, Constants)

        player.payoff_red_zero = round(payoff_red_zero)
        player.payoff_blue_zero = round(payoff_blue_zero)
        player.payoff_red_half = round(payoff_red_half)
        player.payoff_blue_half = round(payoff_blue_half)

        blue_neighbors_half = degree // 2
        red_neighbors_half = degree - blue_neighbors_half

        tries_left = max(Constants.max_retries - player.comprehension_retries, 0)

        return dict(
            role=role,
            degree=degree,
            table_data=table_data,
            blue_neighbors_half=blue_neighbors_half,
            red_neighbors_half=red_neighbors_half,
            tries_left=tries_left,
        )

    def error_message(player, values):
        # skip error messages if retries exceeded
        if player.comprehension_retries >= Constants.max_retries:
            return

        correct_answers = {
            'q_red_zero': player.payoff_red_zero,
            'q_blue_zero': player.payoff_blue_zero,
            'q_red_half': player.payoff_red_half,
            'q_blue_half': player.payoff_blue_half,
        }

        labels = {
            'q_red_zero': '<b>A</b>',
            'q_blue_zero': '<b>B</b>',
            'q_red_half': '<b>C</b>',
            'q_blue_half': '<b>D</b>',
        }

        incorrect_fields = [
            labels[f] for f, v in correct_answers.items()
            if values.get(f) != v
        ]

        if incorrect_fields:
            player.comprehension_retries += 1
            role = player.participant.role

            explanation = (
                "Your payoff only depends on your own shirt choice (Table 1)."
                if role == Constants.minority
                else
                "Your total points = reward for picking a color "
                "(Table 1) + reward for matching neighbors (Table 2)."
            )

            tries_left = max(Constants.max_retries - player.comprehension_retries, 0) + 1

            return (
                f"Incorrect answers: {', '.join(incorrect_fields)}. {explanation} "
                f"<b>You have {tries_left} "
                f"{'try' if tries_left == 1 else 'tries'} left</b>."
            )

    @staticmethod
    def before_next_page(player, timeout_happened):
        if player.comprehension_retries >= Constants.max_retries or timeout_happened:
            player.participant.failed_checks = True


page_sequence = [
    IntroductionPage,
    ComprehensionPage, #only turn off for testing purposes.
                 ]