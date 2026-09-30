from otree.api import *


doc = """
After completing the experiment, participants respond to a set of survey questions.
"""


class Constants(BaseConstants):
    name_in_url = 'questionnaire'
    players_per_group = None
    num_rounds = 1
    majority_role = 'Red'


class Subsession(BaseSubsession):
    pass


class Group(BaseGroup):
    pass

class Player(BasePlayer):
    enjoyment = models.IntegerField(
        choices=[[1, 'Not at all'], [2, 'A little'], [3, 'Somewhat'], [4, 'Very'], [5, 'Extremely']],
        label="Did you enjoy the experiment?",
        widget=widgets.RadioSelect
    )

    clarity = models.IntegerField(
        choices=[[1, 'Very unclear'], [2, 'Unclear'], [3, 'Neutral'], [4, 'Clear'], [5, 'Very clear']],
        label="How clear were the instructions?",
        widget=widgets.RadioSelect
    )

    majority = models.IntegerField(
        choices=[
            [1, 'Red'],
            [2, 'Blue'],
            [3, "Don't know"]
        ],
        label= "Earlier, we told you that players either preferred red or blue, "
              "with one group being larger (the majority) and the other smaller (the minority). " 
                "Which group do you think was the majority?",
        widget=widgets.RadioSelect
    )

    strategy = models.LongStringField(
        blank=True,
        label="What strategy did you use when deciding which color to choose? Please describe your reasoning."
    )

    perceived_rq = models.LongStringField(
        blank=True,
        label="What do you think this experiment was about? What question do you think the researchers are trying to answer?"
    )

    comments = models.LongStringField(
        blank=True,
        label="If you have any other comments or feedback, you can write that here."
    )

# PAGES
class Questionnaire(Page):
    form_model = 'player'

    def is_displayed(player):
        return (
            not player.participant.vars.get("exit_early", False)
            and not player.participant.is_dropout
    )

    @staticmethod
    def get_form_fields(player: Player):
        fields = ['enjoyment', 'clarity', 'majority', 'perceived_rq', 'comments']
        if player.participant.role == Constants.majority_role:
            fields.insert(2, 'strategy') # ask question about strategy only to majority members.
        return fields

page_sequence = [Questionnaire]

