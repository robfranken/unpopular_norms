from otree.api import *
import json
import os
import random
from .functions import compute_utility, payoff_table

from settings import (
    title as TITLE,
    majority_role as MAJORITY,
    minority_role as MINORITY,
    s as S,
    e as E,
    z as Z,
    w as W,
    lambda1 as L1,
    lambda2 as L2,
    base_payment as base,
    max_payment as maxp,
    points_per_euro_majority as PPE1,
    points_per_euro_minority as PPE2,
    num_rounds as nrounds,
    p_minority as p_minority,
    testing as TEST,
)

doc = """
Players enter the group formation page, where they wait until sufficient players (of both roles) arrive.
A group is then formed and the game begins.
Those who arrive too late return to Prolific (two routes, depending on whether they were invited or not).
After timing out, players can return to the game. If they time out twice without returning in time, they exit.
"""

class Constants(BaseConstants):
    title = TITLE
    name_in_url = "fashion_dilemma"
    players_per_group = None
    num_rounds = nrounds
    majority = MAJORITY
    minority = MINORITY
    s = S
    e = E
    z = Z
    w = W
    lambda1 = L1
    lambda2 = L2
    introduction_timeout_seconds = 90
    decision_pages_timeout_seconds = 60
    other_pages_timeout_seconds = 15
    points_per_euro_majority = PPE1
    points_per_euro_minority = PPE2
    base_payment = base
    max_payment = maxp

class Subsession(BaseSubsession):
    def creating_session(self):
        if self.round_number == 1:
            for p in self.get_players():
                p.participant.is_dropout = False

            net_condition = self.session.config.get("network_condition")

            if net_condition and "net_spec" not in self.session.vars:
                file_path = os.path.join("networks", f"network_{net_condition}.json")
                with open(file_path, "r") as f:
                    net = json.load(f)
                self.session.vars["net_spec"] = net

            self.session.vars["group_formed"] = False

        else:
            self.group_like_round(1)

class Player(BasePlayer):
    choice = models.BooleanField(
        verbose_name="Make your choice: Will you wear a Blue or a Red T-shirt today?",
        widget=widgets.RadioSelect,
        choices=[(True, "Blue"), (False, "Red")],
    )
    prolific_id = models.StringField(default=str(" "))
    is_dropout = models.BooleanField(initial=False)
    return_choice = models.BooleanField(initial=False)
    bonus = models.FloatField(initial=0) #field to store bonus (points converted to money, minus base pay)
    arrived_waitpage = models.BooleanField(initial=False)
    arrived_grouppage = models.BooleanField(initial=False)
    checked_neighbors = models.BooleanField(initial=False)

class Group(BaseGroup):
    def set_first_stage_earnings(self):
        players = self.get_players()
        for player in players:
            if player.participant.vars.get("exit_early", False):
                player.payoff = 0
                continue

            my_choice = player.choice
            my_node = player.participant.node
            adj_matrix = player.session.vars["net_spec"]["adj_matrix"]

            neighbors = []
            for i, connection in enumerate(adj_matrix[my_node]):
                if connection == 1:
                    neighbor_player = next(
                        p for p in players if p.participant.node == i
                    )
                    if not neighbor_player.participant.vars.get("exit_early", False):
                        neighbors.append(i)

            neighbor_choices = []
            for neighbor_id in neighbors:
                neighbor_player = next(
                    p for p in players if p.participant.node == neighbor_id
                )
                neighbor_choices.append(neighbor_player.choice)

            utility = compute_utility(
                player_choice=my_choice,
                player_role=player.participant.role,
                neighbors_choices=neighbor_choices,
            )

            player.payoff = max(utility, 0)

def timeout_check(player, timeout_happened):
    participant = player.participant
    if timeout_happened and not participant.is_dropout:
        participant.is_dropout = True
        player.is_dropout = True

def timeout_time(player, timeout_seconds):
    participant = player.participant
    if participant.vars.get('exit', False):
        return 1
    else:
        return timeout_seconds

def group_by_arrival_time_method(subsession, waiting_players):
    session = subsession.session
    group_size = session.config["group_size"]

    if "net_spec" not in session.vars:
        net_condition = session.config.get("network_condition")
        if net_condition:
            file_path = os.path.join("networks", f"network_{net_condition}.json")
            try:
                with open(file_path, "r") as f:
                    net = json.load(f)
                session.vars["net_spec"] = net
            except Exception as e:
                return

    if session.vars.get("group_formed", False):
        for p in waiting_players:
            p.participant.vars["exit_early"] = True
            p.participant.is_dropout = True
        return waiting_players

    net_spec = session.vars.get("net_spec", None)
    adj_matrix = net_spec["adj_matrix"]
    role_vector = net_spec["role_vector"]

    n = len(role_vector)

    role_for_idx = [
        Constants.minority if v == 1 else Constants.majority
        for v in role_vector
    ]

    by_role = {
        Constants.majority: [
            p for p in waiting_players if p.participant.role == Constants.majority
        ],
        Constants.minority: [
            p for p in waiting_players if p.participant.role == Constants.minority
        ],
    }

    required_counts = {
        Constants.majority: sum(1 for r in role_for_idx if r == Constants.majority),
        Constants.minority: sum(1 for r in role_for_idx if r == Constants.minority),
    }

    have_counts = {
        Constants.majority: len(by_role[Constants.majority]),
        Constants.minority: len(by_role[Constants.minority]),
    }
    if (
        have_counts[Constants.majority] >= required_counts[Constants.majority]
        and have_counts[Constants.minority] >= required_counts[Constants.minority]
    ):
        players_ordered = []
        buckets = {
            Constants.majority: by_role[Constants.majority][:],
            Constants.minority: by_role[Constants.minority][:],
        }

        def assign_nodes_and_matrix(selected_players, adj_matrix):

            for i, p in enumerate(selected_players):
                p.participant.node = i
                p.participant.is_dropout = False

        for i in range(n):
            needed_role = role_for_idx[i]
            chosen_player = buckets[needed_role].pop(0)
            players_ordered.append(chosen_player)

        assign_nodes_and_matrix(players_ordered, adj_matrix)

        session.vars["group_formed"] = True
        return players_ordered

    else:
        return

class NetworkFormationWaitPage(WaitPage):
    template_name = "unpop/GroupFormationPage.html"
    group_by_arrival_time = True

    @staticmethod
    def is_displayed(player):
        return (
                player.round_number == 1
                and not player.participant.vars.get("dropout", False)
        )

    def vars_for_template(player):
        if not player.arrived_grouppage:
            player.arrived_grouppage = True

        waiting_players = player.subsession.get_players()
        total_arrived = sum(p.arrived_grouppage for p in waiting_players)

        group_size = player.session.config.get("group_size", len(waiting_players))
        total_needed = int(group_size * 1.3)

        if total_needed == 0:
            percent = 0
            return dict(percent=percent)

        percent = (total_arrived / total_needed) * 100
        percent = min(int(percent), 99)
        return dict(percent=percent, role = player.participant.role)

    @staticmethod
    def after_all_players_arrive(group):
        pass

class IntroductionPage(Page):
    def vars_for_template(player):
        adj_matrix = player.session.vars["net_spec"]["adj_matrix"]
        my_node = player.participant.node
        degree = sum(adj_matrix[my_node])
        table_data = payoff_table(degree)
        group_size = player.session.config["group_size"]

        return dict(
            role=player.participant.role,
            network_condition=player.session.config.get("network_condition"),
            group_size=group_size,
            others=group_size-1,
            degree=degree,
            range_neighbors=list(range(degree + 1)) if degree > 0 else [],
            table_data=table_data,
            base="{:.2f}".format(Constants.base_payment),
            max="{:.2f}".format(Constants.max_payment),
        )

    def is_displayed(player):
        return (
            player.round_number == 1
            and not player.participant.vars.get("exit_early", False)
            and not player.participant.is_dropout
        )

    def get_timeout_seconds(player):
        return timeout_time(player, Constants.introduction_timeout_seconds)

    def before_next_page(player, timeout_happened):
        timeout_check(player, timeout_happened)

class DecisionPage(Page):
    form_model = "player"
    form_fields = ["choice", "checked_neighbors"]

    def get_timeout_seconds(player):
        if player.round_number <= 2:
            timeout = Constants.decision_pages_timeout_seconds
        else:
            timeout = Constants.decision_pages_timeout_seconds / 3
        return timeout_time(player, timeout)

    def before_next_page(player, timeout_happened):
        if 'consecutive_timeouts' not in player.participant.vars:
            player.participant.vars['consecutive_timeouts'] = 0

        if 'exit' not in player.participant.vars:
            player.participant.vars['exit'] = False

        timeout_check(player, timeout_happened)

        if timeout_happened:
            player.participant.vars['consecutive_timeouts'] += 1
        else:
            player.participant.vars['consecutive_timeouts'] = 0
            player.participant.is_dropout = False

        #if playres miss 2 decision rounds and don't return, they exit
        if player.participant.vars['consecutive_timeouts'] >= 2:
            player.participant.vars['exit'] = True

        # auto-submit on time-out
        if timeout_happened or player.participant.is_dropout:
            if player.participant.role == Constants.minority:
                player.choice = True  # always True
            else:
                if player.round_number > 1:
                    prev_choice = player.in_round(player.round_number - 1).choice
                    if prev_choice is not None:  # try to copy previous choice
                        player.choice = prev_choice
                    else:
                        player.choice = (random.random() < p_minority)
                else:  # first round, use probability
                    player.choice = (random.random() < p_minority)

    @staticmethod
    def is_displayed(player: Player):
        return not player.participant.vars.get("exit_early", False)

    def vars_for_template(player):
        adj_matrix = player.session.vars["net_spec"]["adj_matrix"]
        my_node = player.participant.node
        degree = sum(adj_matrix[my_node])

        table_data = payoff_table(degree)

        num_blue_previous_round = 0
        num_red_previous_round = 0
        if player.round_number > 1:
            neighbors = [
                i for i, connection in enumerate(adj_matrix[my_node]) if connection == 1
            ]

            prev_round = player.round_number - 1
            num_blue_previous_round = sum(
                1
                for p in player.group.get_players()
                if p.participant.node in neighbors
                and not p.participant.vars.get("exit_early", False)
                and p.in_round(prev_round).choice is True
            )
            num_red_previous_round = sum(
                1
                for p in player.group.get_players()
                if p.participant.node in neighbors
                and not p.participant.vars.get("exit_early", False)
                and p.in_round(prev_round).choice is False
            )

        return dict(
            group_size=player.session.config["group_size"],
            network_condition=player.session.config.get("network_condition"),
            role=player.participant.role,
            round_number=player.round_number,
            degree=degree,
            range_neighbors=list(range(degree + 1)),
            table_data=table_data,
            num_blue_previous_round=num_blue_previous_round,
            num_red_previous_round=num_red_previous_round,
            is_drop_out=player.participant.is_dropout,
            exit= player.participant.vars.get("exit"),
            timeout_counter = player.participant.vars.get('consecutive_timeouts', 0),
            timeout_seconds = int(DecisionPage.get_timeout_seconds(player))
        )

class ReturnPage(Page):
    form_model = 'player'
    form_fields = ['return_choice']

    def get_timeout_seconds(player):
        return 10

    def is_displayed(player):
        return (
                player.participant.is_dropout
                and not player.participant.vars.get("exit_early", False)
                and not player.participant.vars.get("exit", False)
        )

    def before_next_page(player, timeout_happened):
        if player.return_choice:
            player.participant.is_dropout = False
            # reset consecutive timeouts on return
            player.participant.vars['consecutive_timeouts'] = 0

    def vars_for_template(player):
        return dict(
            exit = player.participant.vars.get("exit"),
            timeout_counter = player.participant.vars.get('consecutive_timeouts', 0),
        )

class ResultsWaitPage(WaitPage):
    template_name = "unpop/ResultsWaitPage.html"

    @staticmethod
    def is_displayed(player: Player):
        return  not player.participant.vars.get("exit_early", False)

    def vars_for_template(player):
        if not player.arrived_waitpage:
            player.arrived_waitpage = True

        players = player.group.get_players()
        arrived = sum(p.arrived_waitpage for p in players)
        total = len(players)

        percent = 100 * arrived / total if total > 0 else 0

        return dict(
            arrived=arrived,
            total=total,
            percent=percent,
            is_drop_out=player.participant.is_dropout,
            exit=player.participant.vars.get("exit"),
        )

    def after_all_players_arrive(group):
        group.set_first_stage_earnings()

    @staticmethod
    def get_timeout_seconds(player):
        return timeout_time(player, 5)

class ResultsPage(Page):
    def vars_for_template(player):
        my_choice = player.choice
        my_payoff = player.payoff

        my_node = player.participant.node
        adj_matrix = player.session.vars["net_spec"]["adj_matrix"]

        neighbors = []
        for i, connection in enumerate(adj_matrix[my_node]):
            if connection == 1:
                neighbors.append(i)

        neighbors_info = []
        for idx, neighbor_id in enumerate(neighbors, start=1):
            neighbor_player = next(
                (
                    p
                    for p in player.group.get_players()
                    if p.participant.node == neighbor_id
                ),
                None,
            )
            if neighbor_player:
                if neighbor_player.choice is None:
                    choice_display = "Missing"
                else:
                    choice_display = "Blue" if neighbor_player.choice else "Red"

                neighbors_info.append(
                    {
                        "neighbor": idx,
                        "id": neighbor_player.id,
                        "choice": choice_display,
                        "payoff": neighbor_player.payoff,
                    }
                )

        my_choice_display = "Blue" if my_choice else "Red"

        return dict(
            my_choice=my_choice_display,
            my_payoff=my_payoff,
            neighbors_info=neighbors_info,
            role=player.participant.role,
            round_number=player.round_number,
        )

    def is_displayed(player):
        return not player.participant.vars.get("exit_early", False) and not player.participant.is_dropout

    def get_timeout_seconds(player):
        return timeout_time(player, Constants.other_pages_timeout_seconds)

class FinalGameResults(Page):
    @staticmethod
    def is_displayed(player):
        return (
            player.round_number == Constants.num_rounds
            and not player.participant.vars.get("exit_early", False)
            and not player.participant.is_dropout
        )

    @staticmethod
    def js_vars(player):
        return dict(completionlink=player.subsession.session.config["completionlink"])

    @staticmethod
    def vars_for_template(player):
        accumulated_earnings = player.participant.payoff
        base = Constants.base_payment

        conversion = (
            Constants.points_per_euro_majority
            if player.participant.role == Constants.majority
            else Constants.points_per_euro_minority
        )

        euros = float(accumulated_earnings) / conversion
        euros = min(euros, Constants.max_payment)
        euros = max(euros, Constants.base_payment)
        bonus = max(euros - base, 0)

        player.participant.vars['bonus'] = round(bonus, 2)

        return dict(
            accumulated_earnings=accumulated_earnings,
            raw_euros=float(accumulated_earnings) / conversion,
            base="{:.2f}".format(base),
            bonus="{:.2f}".format(bonus),
            euros="{:.2f}".format(euros),
            test=TEST,
        )

    @staticmethod
    def before_next_page(player, timeout_happened):
        player.bonus = player.participant.vars.get('bonus', 0)

class ExitPage(Page):
    @staticmethod
    def is_displayed(player: Player):
        return player.participant.vars.get("exit_early", False)

    @staticmethod
    def vars_for_template(player: Player):
        invited = player.participant.vars.get("invited", False)

        if invited:
            message = (
                "Unfortunately, the group for this session is already full. "
                "You will not be participating in the experiment. "
                "You will still receive the base payment for your time and effort."
            )
        else:
            message = (
                "Unfortunately, the group for this session is already full. "
                "You will not be participating in the experiment. "
                "Please return your submission."
            )

        return dict(message=message)

    @staticmethod
    def js_vars(player: Player):
        invited = player.participant.vars.get("invited", False)

        if invited:
            completionlink = player.subsession.session.config.get("completionlink_full")
        else:
            completionlink = player.subsession.session.config.get("completionlink_no_invite")

        return dict(completionlink=completionlink)

page_sequence = [
    NetworkFormationWaitPage,
    IntroductionPage,
    DecisionPage,
    ReturnPage,
    ResultsWaitPage,
    ResultsPage,
    FinalGameResults,
    ExitPage,
]