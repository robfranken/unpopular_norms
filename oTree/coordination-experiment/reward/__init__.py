from otree.api import *


doc = """
Finally, participants are redirected to Prolific with a completion code.
"""


class Constants(BaseConstants):
    name_in_url = 'reward'
    players_per_group = None
    num_rounds = 1

class Subsession(BaseSubsession):
    pass

class Group(BaseGroup):
    pass

class Player(BasePlayer):
    bonus = models.FloatField(initial=0)

# PAGES
class PaymentInfo(Page):
    form_model = 'player'

    @staticmethod
    def js_vars(player):
        if player.participant.is_dropout:
            completionlink = player.subsession.session.config['completionlink_full'] # those who timed out get just the base pay, just like those who couldnt be grouped
        else:
            completionlink = player.subsession.session.config['completionlink']

        return dict(completionlink=completionlink)

page_sequence = [PaymentInfo]
