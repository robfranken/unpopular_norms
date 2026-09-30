from otree.api import *


doc = """
Return to Prolific
"""


class C(BaseConstants):
    NAME_IN_URL = 'return'
    PLAYERS_PER_GROUP = None
    NUM_ROUNDS = 1


class Subsession(BaseSubsession):
    pass


class Group(BaseGroup):
    pass

class Player(BasePlayer):
    confirm = models.BooleanField(
        label="I confirm the above information",
        blank=False
    )


# PAGES
class Info(Page):
    form_model = 'player'
    form_fields = ['confirm']

    @staticmethod
    def is_displayed(player: Player):
        return (
            player.participant.participate is True
            and player.participant.vars.get('consent') is True
            and player.participant.vars.get('failed_checks') is not True
        )

    @staticmethod
    def vars_for_template(player: Player):
            return dict(
                timeslot=player.participant.time_slots,
                participate=player.participant.participate,
                role=player.participant.vars.get('role'),
                #failed_checks=player.participant.vars.get('failed_checks')
            )

class Return(Page):
    @staticmethod
    def js_vars(player):
        session_config = player.subsession.session.config
        if player.field_maybe_none("confirm"):
            link = session_config['participate_completion']
        else:
            link = session_config['no_participate_completion']

        failed = player.participant.vars.get('failed_checks', False)

        return dict(
            completionlink=link,
            failed_checks=bool(failed)
        )

page_sequence = [Info, Return]
