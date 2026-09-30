from otree.api import *


doc = """
Introduction to the experiment
"""


class C(BaseConstants):
    NAME_IN_URL = 'intro'
    PLAYERS_PER_GROUP = None
    NUM_ROUNDS = 1


class Subsession(BaseSubsession):
    pass


class Group(BaseGroup):
    pass


class Player(BasePlayer):
    participate = models.BooleanField(
        label="Do you want to participate in the study?",
        choices=[[True, "Yes"], [False, "No"]],
        widget=widgets.RadioSelect
    )

    time_slots = models.StringField(
        choices=[
            "14:00–14:30 Europe/Amsterdam (13:00–13:30 UK time)"
            #"13:30–14:30",
            #"15:00–16:00",
        ],
        label="Below you find the available time slot(s) for participation.\n"
              "Please only register if you are confident that you can attend the session at the selected time.\n"
              "You will be invited for the session via Prolific.",
        blank=True
    )

    opt_out_reason = models.StringField(
        choices=[
            "Time slot not convenient",
            "Study topic not interesting",
            "Study seems too long",
            "Privacy concerns",
            "Other",
        ],
        widget=widgets.RadioSelect,

        label="Could you indicate your (main) reason?"
    )

    opt_out_reason_text = models.LongStringField(
        label="Please specify:",
        blank=True
    )

    prolific_id = models.StringField(default=str(" "))

# PAGES
class Introduction(Page):
    def before_next_page(player, timeout_happened):
        player.prolific_id = player.participant.label


class Register(Page):
    form_model = 'player'
    form_fields = ['participate', 'time_slots']

    @staticmethod
    def error_message(player, values):
        if values['participate'] and not values['time_slots']:
            return "Please select a time slot if you want to participate."

    @staticmethod
    def before_next_page(player, timeout_happened):
        participant = player.participant
        participant.participate = player.participate
        participant.time_slots = player.time_slots

class OptOut(Page):
    form_model = 'player'
    form_fields = ['opt_out_reason', 'opt_out_reason_text']

    @staticmethod
    def is_displayed(player):
        return player.participate is False



page_sequence = [Introduction, Register, OptOut]
