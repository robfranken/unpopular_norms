import datetime, random
from otree.api import *



doc = """
Consent form for those who want to participate in the experiment
"""

class C(BaseConstants):
    NAME_IN_URL = 'consent'
    PLAYERS_PER_GROUP = None
    NUM_ROUNDS = 1
    BLUE_SLOTS = 10 # the number of participants that are assigned blue (after consenting); oversample to circumvent bottlenecks due to no-show.


class Subsession(BaseSubsession):
    pass


class Group(BaseGroup):
    pass


class Player(BasePlayer):
    consent = models.BooleanField(
        label="",
        choices=[[True, "I consent"], [False, "I don't consent"]],
        blank=False
    )
    consent_timestamp = models.StringField(blank=True)


class ConsentPage(Page):
    form_model = 'player'
    form_fields = ['consent']

    @staticmethod
    def is_displayed(player: Player):
        return player.participant.participate is True

    @staticmethod
    def before_next_page(player: Player, timeout_happened):
        if player.consent:
            player.participant.vars['consent'] = True
            player.consent_timestamp = datetime.datetime.now().isoformat()

            # assign role from quota:
            session = player.session
            blue_count = session.vars.get("blue_count", 0)
            if blue_count < C.BLUE_SLOTS:
                role = "Blue"
                session.vars["blue_count"] = blue_count + 1
            else:
                role = "Red"
            player.participant.role = role

page_sequence = [ConsentPage]
