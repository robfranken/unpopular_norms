from otree.api import Bot
from . import *


class PlayerBot(Bot):
    def play_round(self):
        # Introduction Page
        if self.round_number == 1:
            yield IntroductionPage

        # comprehension check
        player = self.player
        #if self.round_number == 1 and player.participant.role == Constants.majority_role:
        #    yield ComprehensionPage, {
        #        'q_red_zero': player.payoff_red_zero,
        #        'q_blue_zero': player.payoff_blue_zero,
        #        'q_red_half': player.payoff_red_half,
        #        'q_blue_half': player.payoff_blue_half,
        #    }

        # DecisionPage: players follow their preference
        choice = True if self.player.participant.role == Constants.minority_role else False
        yield DecisionPage, {'choice': choice}

        # results
        yield ResultsPage

        # Final page in last round
        if self.round_number == Constants.num_rounds:
            if not self.player.group.failed and not self.player.participant.is_dropout:
                yield FinalGameResults
            else:
                yield FailedGamePage


