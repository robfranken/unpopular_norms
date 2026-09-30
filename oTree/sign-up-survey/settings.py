from os import environ

SESSION_CONFIGS = [
     dict(
         name='sign_up',
         app_sequence=['intro', 'consent', 'comprehension', 'return'],
         num_demo_participants=200,
         participate_completion='https://app.prolific.com/submissions/complete?cc=C10D93TZ',
         no_participate_completion='https://app.prolific.com/submissions/complete?cc=C1NPREPM',
     ),
]

SESSION_CONFIG_DEFAULTS = dict(
    real_world_currency_per_point=1.00, participation_fee=0.00, doc=""
)

PARTICIPANT_FIELDS = ["participate", "time_slots", "consent", "role", "failed_checks"]
SESSION_FIELDS = []

ROOMS = [
    dict(
        name='signup',
        display_name='Sign-up survey for interactive group experiment',
    ),
]


LANGUAGE_CODE = 'en'
REAL_WORLD_CURRENCY_CODE = 'USD'
USE_POINTS = True

ADMIN_USERNAME = 'admin'
ADMIN_PASSWORD = environ.get('OTREE_ADMIN_PASSWORD')

DEMO_PAGE_INTRO_HTML = """ """

SECRET_KEY = '9795047743176'
