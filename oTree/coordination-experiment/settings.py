from os import environ

SESSION_CONFIGS = [
    dict(
        name="unpopular_norm_4",
        display_name="test_n4",
        num_demo_participants=200,
        group_size=4,
        network_condition="test_n4",
        app_sequence=[ "arrive", "unpop", "survey", "reward"],
        completionlink='https://app.prolific.com/submissions/complete?cc=CI5BFLAB',
        completionlink_full='https://app.prolific.com/submissions/complete?cc=CKGTTDJJ',
        completionlink_no_invite='https://app.prolific.com/submissions/complete?cc=CG1SAGQG',
        use_browser_bots=False,
    ),

    dict(
        name="unpopular_norm_50_control",
        display_name="test_n50_control",
        num_demo_participants=200,
        group_size=50,
        network_condition="test_n50_random",
        app_sequence=[ "arrive", "unpop", "survey", "reward"],
        completionlink='https://app.prolific.com/submissions/complete?cc=CI5BFLAB',
        completionlink_full='https://app.prolific.com/submissions/complete?cc=CKGTTDJJ',
        completionlink_no_invite='https://app.prolific.com/submissions/complete?cc=CG1SAGQG',
        use_browser_bots=False,
    ),


    dict(
        name="unpopular_norm_50",
        display_name="test_n50",
        num_demo_participants=200,
        group_size=50,
        network_condition="test_n50",
        app_sequence=[ "arrive", "unpop", "survey", "reward"],
        completionlink='https://app.prolific.com/submissions/complete?cc=CI5BFLAB',
        completionlink_full='https://app.prolific.com/submissions/complete?cc=CKGTTDJJ',
        completionlink_no_invite='https://app.prolific.com/submissions/complete?cc=CG1SAGQG',
        use_browser_bots=False,
    ),
]

# set some central parameters to be used across apps:
title = 'The Fashion Dilemma'
majority_role = 'Red'
minority_role = 'Blue'
p_minority = 0.1 # !!this needs to correspond to the proportion of minorities in the network configuration!!
num_rounds = 25

# including also the incentive structure
s = 15
e = 10
z = 50
w = 40
lambda1 = 4.3
lambda2 = 1.8

# and payment variables (base pay; conversion rates; etc.)
base_payment = 2.5 #base pay of 2.50 (for estimated 25 min.)
max_payment = 7.5 #max of 7.50
points_per_euro_majority = 167  #200 (assuming 30 rounds)
points_per_euro_minority = 33  #40 (assuming 30 rounds)

#configure a room
ROOMS = [
    dict(
        name='1',
        display_name='Network #1: heterogenous (central fanatics)',
        #participant_label_file='_rooms/fashion_dilemma.txt',
        #use_secure_urls=True,
    ),
    dict(
        name='2',
        display_name='Network #2: homogenous (random)',
        #participant_label_file='_rooms/fashion_dilemma.txt',
        #use_secure_urls=True,
    ),
]

SESSION_CONFIG_DEFAULTS = dict(
    real_world_currency_per_point=1/30,
    participation_fee=3.00,
    doc="",
)

PARTICIPANT_FIELDS = [ "bonus", "consent", "is_dropout", "consecutive_timeouts", "exit", "role", 'has_dropped_out', 'node', 'adj_matrix', 'role_vector', 'exit_early', 'failed_checks']
LANGUAGE_CODE = "en"
REAL_WORLD_CURRENCY_CODE = "EUR"
USE_POINTS = True
ADMIN_USERNAME = "admin"
ADMIN_PASSWORD = environ.get("OTREE_ADMIN_PASSWORD")
DEMO_PAGE_INTRO_HTML = """ """
SECRET_KEY = "secret"