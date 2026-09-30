from . import *

class PlayerBot(Bot):
    def play_round(self):
        answers = {
            'enjoyment': 4,  # "Very"
            'clarity': 4,    # "Clear"
            'majority': 1,
            'perceived_rq': "I think the experiment studies social coordination.",
            'comments': "No additional comments."
        }

        if self.player.participant.role == Constants.majority_role:
            answers['strategy'] = "I tried to match what I thought most others would choose."

        yield Questionnaire, answers
