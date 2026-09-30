from otree.api import *
import csv
import os
import random

doc = """
Invited players arrive,
their Prolific ID is read
And the linked role is retrieved from the key file participants.csv
"""

# key file containing (invited) Prolific IDs and assigned roles
# from 'sign-up-experiment'
PARTICIPANT_ROLES = {}
csv_path = os.path.join(os.path.dirname(__file__), 'participants.csv')

with open(csv_path, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        label = row["participant.label"].strip()
        role = row["participant.role"].strip()
        PARTICIPANT_ROLES[label] = role

class C(BaseConstants):
    NAME_IN_URL = 'arrive'
    PLAYERS_PER_GROUP = None
    NUM_ROUNDS = 1

class Subsession(BaseSubsession):
    pass

class Group(BaseGroup):
    pass

class Player(BasePlayer):
    prolific_id = models.StringField(default="")

class Introduction(Page):
    def before_next_page(player, timeout_happened):
        # read in Prolific ID
        participant_label = player.participant.label
        player.prolific_id = participant_label

        # get assigned role (and if none assigned--that is, non-inivited player--give random)
        role = PARTICIPANT_ROLES.get(participant_label)

        if role:
            player.participant.vars['role'] = role
            player.participant.vars['invited'] = True

        else:
            session = player.session
            role_pool = ["Red", "Blue"]
            if "role_index" not in session.vars:
                session.vars["role_index"] = 0
            idx = session.vars["role_index"]
            assigned_role = role_pool[idx % len(role_pool)]
            session.vars["role_index"] += 1

            player.participant.vars['role'] = assigned_role
            player.participant.vars['invited'] = False # non-invited participants don't receive a show-up fee...

page_sequence = [Introduction]