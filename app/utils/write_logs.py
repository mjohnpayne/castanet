import json
import os


def write_input_params(payload):
    try:
        payload_for_sav = payload.copy()
        for key, val in payload_for_sav.items():
            '''Frozensets don't serialise to JSON'''
            if type(val) == frozenset:
                payload_for_sav[key] = list(val)
        '''Make the directories if they don't exist (they really shouldn't, but y'kna, reruns)'''
        if not os.path.exists(f"{payload['SaveDir']}/"):
            os.mkdir(f"{payload['SaveDir']}/")
        if not os.path.exists(f"{payload['SaveDir']}/{payload['ExpName']}/"):
            os.mkdir(f"{payload['SaveDir']}/{payload['ExpName']}/")
        with open(f"{payload['SaveDir']}/{payload['ExpName']}/run_parameters.json", "w") as f:
            f.write(json.dumps(payload_for_sav))
    except FileNotFoundError as e:
        print("INFO: Couldn't handle logging for some reason. This really should work, please let us know if this borks.")
